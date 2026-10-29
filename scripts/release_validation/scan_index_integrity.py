#!/usr/bin/env python3
"""Sequential full-file SHA256 scans. Never change an input; continue after a failed file.

A failure offset is the exact start/length of the failing read, not a claim
that the first byte in that range is the first damaged byte on disk.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import time
import subprocess
import sys

from verify_inputs import SUFFIXES, metadata


def scan_file(path, chunk_size=4*1024*1024):
    path = Path(path)
    start = time.monotonic()
    result = dict(path=str(path), status='UNREADABLE', sha256=None)
    h, offset = hashlib.sha256(), 0
    try:
        before = path.stat()
        result['size'] = before.st_size
        with path.open('rb', buffering=0) as f:
            while offset < before.st_size:
                size = min(chunk_size, before.st_size-offset)
                try:
                    data = f.read(size)
                except OSError as e:
                    result.update(failure_offset=offset, requested_bytes=size, errno=e.errno, error=str(e))
                    break
                if not data:
                    result.update(failure_offset=offset, requested_bytes=size, error='unexpected EOF')
                    break
                h.update(data)
                offset += len(data)
            else:
                after = path.stat()
                if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
                    result.update(error='file changed during read', status='CHANGED')
                else:
                    result.update(status='READABLE', sha256=h.hexdigest())
    except OSError as e:
        result.update(failure_offset=offset, errno=e.errno, error=str(e))
    result.update(bytes_completely_read=offset, seconds=round(time.monotonic()-start, 3))
    return result


def bounded_scan(path, timeout=90):
    command = [sys.executable, str(Path(__file__).resolve()), '--single-file', str(path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
        if process.returncode:
            return dict(path=str(path), status='SCAN_ERROR', sha256=None, error=stderr)
        return json.loads(stdout)
    except subprocess.TimeoutExpired:
        offset = None
        # Linux diagnostics: this is the offset of the blocked read, not a
        # claim about the first physically corrupt byte.
        fdroot = Path('/proc')/str(process.pid)/'fd'
        for fd in fdroot.glob('*'):
            try:
                if fd.resolve() == Path(path).resolve():
                    info = fdroot.parent.joinpath('fdinfo', fd.name).read_text()
                    offset = int(next(l.split()[1] for l in info.splitlines() if l.startswith('pos:')))
            except (OSError, StopIteration):
                pass
        process.kill()
        process.communicate()
        return dict(path=str(path), status='TIMEOUT', sha256=None,
                    failure_offset=offset, requested_bytes=4*1024*1024,
                    error='full-file read exceeded {} seconds; not fully verified'.format(timeout))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest')
    p.add_argument('--single-file')
    p.add_argument('--resume-jsonl')
    p.add_argument('--read-timeout', type=int, default=90)
    p.add_argument('--output')
    a = p.parse_args()
    if a.single_file:
        print(json.dumps(scan_file(a.single_file)))
        return
    if not a.manifest or not a.output:
        p.error('--manifest and --output required')
    previous = {}
    if a.resume_jsonl:
        for line in Path(a.resume_jsonl).read_text().splitlines():
            row = json.loads(line)
            previous[row['path']] = row
    out = Path(a.output)
    out.mkdir(exist_ok=False)
    rows = list(csv.DictReader(open(a.manifest), delimiter='\t'))
    summary = dict(io_concurrency=1, chunk_bytes=4*1024*1024, assemblies=[])
    for row in rows:
        asm = row['assembly_id']
        entry = dict(assembly=asm, components=[])
        try:
            entry['metadata'] = metadata(row)
        except Exception as e:
            entry['metadata'] = dict(status='FAIL', error=str(e))
        paths = [row['index_prefix']+'.'+s for s in SUFFIXES]
        # Include exact source FASTA/FAI/GZI so rebuilding can be approved on
        # completely readable, frozen source bytes rather than filenames alone.
        paths += [row['fasta']+s for s in ['', '.fai', '.gzi']]
        for path in paths:
            print(asm+' reading '+path, flush=True)
            if path in previous:
                result = dict(previous[path], reused_from=a.resume_jsonl)
                result.pop('assembly', None)
            else:
                result = bounded_scan(path, a.read_timeout)
            entry['components'].append(result)
            with (out/'files.jsonl').open('a') as h:
                h.write(json.dumps(dict(assembly=asm, **result))+'\n')
            print(asm+' '+result['status']+' '+str(result.get('failure_offset', '')),
                  flush=True)
        entry['status'] = ('READABLE_METADATA_MATCH' if
                           entry['metadata']['status'] == 'METADATA_PASS' and
                           all(x['status'] == 'READABLE' for x in entry['components']) else 'FAIL')
        summary['assemblies'].append(entry)
        (out/(asm+'.json')).write_text(json.dumps(entry, indent=2)+'\n')
        (out/'summary.partial.json').write_text(json.dumps(summary, indent=2)+'\n')
    summary['status'] = 'PASS' if all(x['status']=='READABLE_METADATA_MATCH' for x in summary['assemblies']) else 'FAIL'
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(summary['status'], flush=True)


if __name__ == '__main__':
    main()
