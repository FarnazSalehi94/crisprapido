#!/usr/bin/env python3
"""Guarded bounded old-vs-new performance benchmark.

Runs ONE configured engine/size per Slurm task under the memory/disk supervisor.
No entry point launches a full library or cohort. Old-engine output is diagnostic
only; it is never published as a valid search result.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2)
        f.write('\n')


def available_kib():
    fields = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    return int(fields['MemAvailable'].split()[0])


def require(condition, message):
    if not condition:
        raise ValueError(message)


REFERENCE_SHA256 = 'bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b'
EXPECTED = dict(format='hprcv2-performance-bench-v1', node='octopus02', assembly='HG00735_pat',
                reference_sha256=REFERENCE_SHA256,
                fai_sha256='7f6187771a76e8af51b983e2df4e00970dd72722dea0b9bd5c425ef817e3c26d',
                settings=dict(m=2, b=1, z=1, f=0.75, pam='GG'), sizes=[1, 100, 1000])
LIMITS = dict(stage_timeout_seconds=3600, min_start_available_kib=67108864,
              min_running_available_kib=8388608, min_scratch_bytes=60000000000,
              min_running_disk_bytes=10000000000, file_limit_bytes=10737418240,
              new_rss_kib=41943040, old_rss_kib=115343360)
SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()


def validate(plan):
    for key, value in EXPECTED.items():
        require(plan.get(key) == value, 'plan drift: '+key)
    require(plan['limits'] == LIMITS, 'resource-guard drift')
    require(plan['library']['sha256'] == 'a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971', 'library drift')
    require(plan['old']['binary_sha256'] == 'ee543e772ac6e956e8fc3eb88e6ea52517e5192eff9da255902c3bad6073faed', 'old binary drift')
    require(plan['old']['columba_sha256'] == '7e44974c81dad4943292fa64a71e54073908b4909ce327480464d769322c40d1', 'old columba drift')
    require(set(plan['old']['components']) == set(SUFFIXES), 'old index component drift')
    require(plan['new']['threads'] == 1 and plan['old']['threads'] == 16, 'thread configuration drift')
    for key in ['worker', 'validator', 'supervisor']:
        require(Path(plan[key]).is_absolute() and str(Path(plan[key])) in plan['pins'], 'unpinned tooling: '+key)
    return plan


def verify_pins(plan):
    for path, expected in plan['pins'].items():
        require(digest(path) == expected, 'SHA256 mismatch: '+path)
    binding = json.loads(Path(plan['new']['binding']).read_text())
    require(binding['reference_sha256'] == plan['reference_sha256'], 'native binding/reference mismatch')


def run(command, out, err, metrics, plan, scratch, phase):
    limits = plan['limits']
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    with Path(out).open('x') as o, Path(err).open('x') as e:
        proc = subprocess.Popen(list(map(str, command)), stdout=o, stderr=e, start_new_session=True)
        try:
            while proc.poll() is None:
                require(time.monotonic()-start < limits['stage_timeout_seconds'], phase+' timeout')
                require(available_kib() >= limits['min_running_available_kib'], phase+' available-memory guard')
                require(min(shutil.disk_usage(scratch).free, shutil.disk_usage(plan['root']).free)
                        >= limits['min_running_disk_bytes'], phase+' free-space guard')
                time.sleep(1)
            require(proc.returncode == 0, phase+' returned '+str(proc.returncode))
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    save(metrics, dict(phase=phase, command=list(map(str, command)),
                      wall_seconds=time.monotonic()-start,
                      user_seconds=after.ru_utime-before.ru_utime,
                      system_seconds=after.ru_stime-before.ru_stime,
                      cumulative_max_individual_descendant_rss_kib=after.ru_maxrss,
                      caveat='ru_maxrss is cumulative maximum individual child, not aggregate; tree peak is supervised separately.'))


def parse_native(text):
    found = re.findall(r'^states=(\d+) search_wall_seconds=([0-9.eE+-]+)$', text, re.M)
    require(len(found) == 1, 'missing/duplicate native timing line')
    return dict(states=int(found[0][0]), traversal_seconds=float(found[0][1]))


def parse_time_v(text):
    def grab(pattern):
        m = re.search(pattern, text, re.M)
        require(m is not None, 'time -v field missing: '+pattern)
        return m.group(1)
    wall = grab(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)')
    parts = wall.strip().split(':')
    require(len(parts) in (2, 3) and all(p.replace('.', '', 1).isdigit() for p in parts), 'unparsed wall format')
    seconds = sum(float(p)*60**i for i, p in enumerate(reversed(parts)))
    result = dict(wall_seconds=seconds,
                  cpu_percent=float(grab(r'Percent of CPU this job got: ([0-9.]+)%')),
                  user_seconds=float(grab(r'User time \(seconds\): ([0-9.]+)')),
                  system_seconds=float(grab(r'System time \(seconds\): ([0-9.]+)')),
                  max_rss_kib=int(grab(r'Maximum resident set size \(kbytes\): (\d+)')))
    require(int(grab(r'Exit status: (\d+)')) == 0, 'old engine nonzero exit')
    return result


def copy_checked(source, dest, expected):
    h = hashlib.sha256()
    with open(source, 'rb') as src, open(dest, 'xb') as dst:
        for block in iter(lambda: src.read(4*1024*1024), b''):
            dst.write(block)
            h.update(block)
    require(h.hexdigest() == expected, 'staging checksum mismatch: '+str(source))
    return dest


def stage_executable(source, dest, expected):
    """Staged executables keep the executable bit; a 0600 copy is not runnable."""
    path = Path(copy_checked(source, dest, expected))
    path.chmod(0o700)
    return path


def inner(plan, engine, size, attempt):
    verify_pins(plan)
    limits = plan['limits']
    require(available_kib() >= limits['min_start_available_kib'], 'insufficient available memory at start')
    scratch_root = Path('/scratch')/os.environ['USER']
    scratch_root.mkdir(mode=0o700, exist_ok=True)
    require(shutil.disk_usage(scratch_root).free >= limits['min_scratch_bytes'], 'insufficient scratch')
    scratch = Path(tempfile.mkdtemp(prefix='bench-'+os.environ['SLURM_JOB_ID']+'-', dir=str(scratch_root)))
    save(attempt/'preflight.json', dict(host=socket.gethostname(), job=os.environ['SLURM_JOB_ID'],
          scratch=str(scratch), mem_available_kib=available_kib(),
          scratch_free_bytes=shutil.disk_usage(scratch).free, engine=engine, size=size))
    resource.setrlimit(resource.RLIMIT_FSIZE, (limits['file_limit_bytes'], limits['file_limit_bytes']))
    os.environ.update(TMPDIR=str(scratch), PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
    audit_ref = scratch/'audit-reference.fa'
    audit_ref.symlink_to(plan['reference'])
    shutil.copyfile(plan['fai'], str(audit_ref)+'.fai')
    try:
        if engine == 'new':
            guides = plan['guides'][str(size)]
            require(digest(guides) == plan['pins'][guides], 'guide file drift')
            paf = attempt/'new.paf'
            command = [plan['new']['cli'], '--search-backend', 'pam-anchored', '--pam-search-bin', plan['new']['native'],
                       '--columba-index', plan['new']['index'], '-r', plan['reference'], '-G', guides,
                       '-m', 2, '-b', 1, '-z', 1, '-f', '.75', '-p', 'GG', '-t', 1,
                       '--mismatch-scores', plan['matrices']['mismatch'], '--pam-scores', plan['matrices']['pam']]
            run(command, paf, attempt/'search.stderr', attempt/'search-metrics.json', plan, scratch, 'search')
            native = parse_native((attempt/'search.stderr').read_text())
            audit = attempt/'independent-audit'
            command = ['/usr/bin/python3', Path(plan['release_scripts'])/'audit_emitted.py', '--paf', paf,
                       '--fasta', audit_ref, '--guides', guides, '--matrices', Path(plan['matrices']['mismatch']).parent,
                       '--samtools', plan['new']['samtools'], '--output', audit,
                       '--max-mismatches', 2, '--max-bulges', 1, '--max-bulge-size', 1, '--min-match-fraction', '.75']
            run(command, attempt/'audit.stdout', attempt/'audit.stderr', attempt/'audit-metrics.json', plan, scratch, 'audit')
            summary = json.loads((audit/'summary.json').read_text())
            require(summary['status'] == 'PASS' and summary['errors'] == 0, 'new-engine audit failed')
            save(attempt/'summary.json', dict(engine=engine, size=size, guides=size,
                search_wall_seconds=json.loads((attempt/'search-metrics.json').read_text())['wall_seconds'],
                native_states=native['states'], native_traversal_seconds=native['traversal_seconds'],
                rows_emitted=summary['rows'], audit_status='PASS',
                audit_wall_seconds=json.loads((attempt/'audit-metrics.json').read_text())['wall_seconds'],
                threads=plan['new']['threads'], allocated_cpus=plan['new']['cpus'],
                caveat='Single fresh call, quiet sequential node, uncontrolled caches; not a repeated benchmark.'))
        else:
            require(size in plan['sizes'], 'unknown size')
            new_paf = Path(plan['root'])/'01_new_backend'/('N%d' % size)
            new_paf = next(new_paf.glob('attempt-*/new.paf'))
            require(new_paf.is_file(), 'new-engine baseline for containment missing')
            stage = scratch/'input'
            stage.mkdir()
            index_dir = scratch/'index'
            index_dir.mkdir()
            start = time.monotonic()
            staged_ref = copy_checked(plan['reference'], stage/'HG00735_pat.fa', plan['pins'][plan['reference']])
            shutil.copyfile(plan['fai'], str(staged_ref)+'.fai')
            mm = copy_checked(plan['matrices']['mismatch'], stage/'mismatch_scores.txt', plan['pins'][plan['matrices']['mismatch']])
            pm = copy_checked(plan['matrices']['pam'], stage/'pam_scores.txt', plan['pins'][plan['matrices']['pam']])
            guides = copy_checked(plan['guides'][str(size)], stage/'guides.fa', plan['pins'][plan['guides'][str(size)]])
            binary = stage_executable(plan['old']['binary'], stage/'crisprapido_newCFD', plan['old']['binary_sha256'])
            columba = stage_executable(plan['old']['columba'], stage/'columba', plan['old']['columba_sha256'])
            prefix = index_dir/'HG00735_pat'
            for suffix in SUFFIXES:
                copy_checked(plan['old']['index_prefix']+'.'+suffix, str(prefix)+'.'+suffix, plan['old']['components'][suffix])
            staging_seconds = time.monotonic()-start
            lib = plan['old']['library_path']
            paf = attempt/'old.paf'
            command = [plan['old']['time_bin'], '-v', '-o', attempt/'old.time', 'env', '-u', 'LD_PRELOAD',
                       '-u', 'LD_AUDIT', 'LD_LIBRARY_PATH='+lib, plan['old']['loader'], '--library-path', lib,
                       str(binary), '-r', staged_ref, '-G', guides, '-p', 'GG', '-m', 2, '-b', 1, '-z', 1,
                       '-f', '0.75', '-t', plan['old']['threads'], '--mismatch-scores', mm, '--pam-scores', pm,
                       '--columba-bin', columba, '--columba-index', prefix]
            run(command, paf, attempt/'old.stderr', attempt/'old-metrics.json', plan, scratch, 'old-search')
            timing = parse_time_v((attempt/'old.time').read_text())
            validation = attempt/'old-validation'
            command = ['/usr/bin/python3', plan['validator'], '--paf', paf, '--guides', guides,
                       '--reference', audit_ref, '--matrices', Path(plan['matrices']['mismatch']).parent,
                       '--samtools', plan['new']['samtools'], '--release-scripts', plan['release_scripts'],
                       '--new-paf', new_paf, '--out', validation]
            run(command, attempt/'validate.stdout', attempt/'validate.stderr', attempt/'validate-metrics.json', plan, scratch, 'old-validation')
            checked = json.loads((validation/'summary.json').read_text())
            require(checked['rows'] == checked['accounted'], 'old-row validation accounting mismatch')
            save(attempt/'summary.json', dict(engine=engine, size=size, guides=size,
                staging_seconds=staging_seconds, **timing, rows_emitted=checked['rows'],
                strict_eligible_witnesses=checked['strict_eligible'],
                strict_eligible_contained_in_new=checked['strict_eligible_contained_in_new'],
                strict_eligible_missing_in_new=checked['strict_eligible_missing_in_new'],
                all_rows_contained_in_new=checked['rows_contained_in_new'],
                validation=checked, threads=plan['old']['threads'], allocated_cpus=plan['old']['cpus'],
                stderr_bytes=(attempt/'old.stderr').stat().st_size,
                caveat='Preserved historical binary with legacy PAM/geometry semantics; output is a performance/containment diagnostic, not a valid search result.'))
        verify_pins(plan)
        outputs = {str(p): dict(sha256=digest(p), bytes=p.stat().st_size)
                   for p in attempt.rglob('*') if p.is_file() and p.name not in ['monitor.log', 'monitor-summary.json', 'INCOMPLETE.json']}
        save(attempt/'STAGES_COMPLETE.json', dict(status='STAGES_COMPLETE_AWAITING_SUPERVISOR', outputs=outputs))
    finally:
        shutil.rmtree(scratch)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', required=True, type=Path)
    p.add_argument('--engine', required=True, choices=['new', 'old'])
    p.add_argument('--size', required=True, type=int)
    p.add_argument('--inner', type=Path)
    a = p.parse_args()
    plan = validate(json.loads(a.plan.read_text()))
    require(a.size in plan['sizes'], 'size not in plan')
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURMD_NODENAME'), 'Slurm execution required')
    require(os.environ['SLURMD_NODENAME'] == plan['node'], 'wrong node')
    require(int(os.environ['SLURM_CPUS_PER_TASK']) == plan[a.engine]['cpus'], 'CPU allocation drift')
    if a.inner:
        inner(plan, a.engine, a.size, a.inner)
        return
    stage = '01_new_backend' if a.engine == 'new' else '02_old_pipeline'
    logical = Path(plan['root'])/stage/('N%d' % a.size)
    logical.mkdir(parents=True, exist_ok=True)
    with (logical/'task.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        success = logical/'SUCCESS.json'
        if success.exists():
            result = json.loads(success.read_text())
            require(result['plan_sha256'] == digest(a.plan), 'completed task has a different plan')
            for path, spec in result['outputs'].items():
                require(digest(path) == spec['sha256'], 'completed output changed: '+path)
            print('Already complete; checksums verified; no benchmark rerun.')
            return
        attempt = logical/('attempt-'+os.environ['SLURM_JOB_ID'])
        attempt.mkdir()
        save(attempt/'INCOMPLETE.json', dict(plan=str(a.plan), engine=a.engine, size=a.size))
        try:
            verify_pins(plan)
            limit = plan['limits']['new_rss_kib' if a.engine == 'new' else 'old_rss_kib']
            command = list(map(str, ['/usr/bin/python3', plan['supervisor'], '--limit-kib', limit, '--interval', '.1',
                       '--log', attempt/'monitor.log', '--summary', attempt/'monitor-summary.json', '--',
                       '/usr/bin/python3', plan['worker'], '--plan', a.plan, '--engine', a.engine,
                       '--size', a.size, '--inner', attempt]))
            subprocess.run(command, check=True)
            monitor = json.loads((attempt/'monitor-summary.json').read_text())
            require(monitor['reason'] == 'completed', 'supervisor did not pass: '+str(monitor['reason']))
            data = json.loads((attempt/'STAGES_COMPLETE.json').read_text())
            for path, spec in data['outputs'].items():
                require(digest(path) == spec['sha256'], 'output checksum mismatch')
            for name in ['monitor.log', 'monitor-summary.json', 'STAGES_COMPLETE.json']:
                path = attempt/name
                data['outputs'][str(path)] = dict(sha256=digest(path), bytes=path.stat().st_size)
            save(success, dict(status='COMPLETE', plan_sha256=digest(a.plan), attempt=str(attempt),
                               outputs=data['outputs'], monitor=monitor))
            (attempt/'INCOMPLETE.json').unlink()
        except BaseException as exc:
            save(attempt/'FAILED.json', dict(error=str(exc)))
            raise


if __name__ == '__main__':
    main()
