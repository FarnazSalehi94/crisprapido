#!/usr/bin/env python3
"""Fail-closed resource guard for a bounded recovery wave; never submits jobs.

A shared STOP record blocks subsequent tasks. Running guarded tasks terminate
on STOP/resource failure, leaving all recovery evidence in place. Wrap the
existing process-tree supervisor so it receives TERM and cleans descendants.
"""
import argparse
import ctypes
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import time


def available_bytes():
    fields = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    return int(fields['MemAvailable'].split()[0]) * 1024


def resource_error(stop, scratch, storage, scratch_floor, storage_floor, memory_floor):
    if stop.exists():
        return 'wave STOP exists'
    for path, floor in ((scratch, scratch_floor), (storage, storage_floor)):
        if shutil.disk_usage(path).free < floor:
            return 'free-space floor: ' + str(path)
    if available_bytes() < memory_floor:
        return 'available-memory floor'
    return None


def latch(stop, reason):
    # Append rather than overwrite: retain concurrent failure reports.
    with stop.open('a') as f:
        f.write(json.dumps(dict(reason=reason, host=socket.gethostname(),
            job=os.environ.get('SLURM_JOB_ID'), time=time.time())) + '\n')
        f.flush()
        os.fsync(f.fileno())


def become_subreaper():
    if ctypes.CDLL(None, use_errno=True).prctl(36, 1, 0, 0, 0) != 0:
        raise RuntimeError('cannot become child subreaper')


def descendants():
    parents = {}
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():
            continue
        try:
            fields = dict(line.split(':', 1) for line in (path/'status').read_text().splitlines() if ':' in line)
            parents[int(path.name)] = int(fields['PPid'])
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    members = {os.getpid()}
    while True:
        expanded = members | {pid for pid, parent in parents.items() if parent in members}
        if expanded == members:
            return members - {os.getpid()}
        members = expanded


def cleanup(proc, timeout=15):
    # Subreaping and /proc ancestry cover workers that create separate sessions.
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        pass
    for _ in range(4):
        members = descendants()
        for sig in (signal.SIGSTOP, signal.SIGKILL):
            for pid in members:
                try:
                    os.kill(pid, sig)
                except ProcessLookupError:
                    pass
        time.sleep(.05)
    proc.wait(timeout=5)
    while True:
        try:
            pid, _ = os.waitpid(-1, os.WNOHANG)
            if pid == 0:
                break
        except ChildProcessError:
            break


def fail_closed(stop, reason, proc):
    # A broken/full storage mount must never prevent local process cleanup.
    try:
        latch(stop, reason)
    except OSError as exc:
        print('STOP persistence FAILED: ' + str(exc), file=sys.stderr)
    finally:
        if proc is not None:
            cleanup(proc)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--wave-dir', type=Path, required=True)
    p.add_argument('--scratch', type=Path, default=Path('/scratch'))
    p.add_argument('--storage', type=Path, default=Path('/moosefs/salehi/HPRCv2_CRISPRapido'))
    p.add_argument('--scratch-floor', type=int, default=20_000_000_000)
    p.add_argument('--storage-floor', type=int, default=2_000_000_000_000)
    p.add_argument('--memory-floor', type=int, default=8_589_934_592)
    p.add_argument('--interval', type=float, default=5)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    if not command or not a.wave_dir.is_dir() or a.interval <= 0:
        p.error('existing wave directory, command, and positive interval required')
    stop = a.wave_dir/'STOP'
    proc = None
    received = []
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, lambda sig, frame: received.append(sig))
    reason = None
    try:
        become_subreaper()
        reason = resource_error(stop, a.scratch, a.storage, a.scratch_floor,
                                a.storage_floor, a.memory_floor)
        if reason or received:
            raise RuntimeError(reason or ('signal ' + str(received[0])))
        proc = subprocess.Popen(command, start_new_session=True)
        while proc.poll() is None:
            reason = ('signal ' + str(received[0])) if received else resource_error(
                stop, a.scratch, a.storage, a.scratch_floor, a.storage_floor, a.memory_floor)
            if reason:
                raise RuntimeError(reason)
            time.sleep(a.interval)
        if received:
            raise RuntimeError('signal ' + str(received[0]))
        if proc.returncode != 0:
            raise RuntimeError('child exit ' + str(proc.returncode))
        if descendants():
            raise RuntimeError('orphan descendant after child exit')
    except BaseException as exc:
        reason = str(exc)
        fail_closed(stop, reason, proc)
        raise SystemExit(reason)


if __name__ == '__main__':
    main()
