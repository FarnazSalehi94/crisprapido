#!/usr/bin/env python3
"""Restart-safe full-library Brunello worker (PREPARED; requires explicit plan).

One Slurm task per assembly. Guides are processed in checksummed chunks inside
one job so a killed search loses at most one bounded chunk; completed chunks
resume from checksum-bound checkpoints. Every assembly index/reference binding
is fully verified per run. Search, independent audit and neural annotation are
separate stages with separate completion markers. This file launches nothing by
itself: execution requires a plan whose `authorization` token names the
approving user decision, which this code cannot fabricate.
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

HERE = Path(__file__).resolve().parent
for candidate in (HERE, HERE.parent, HERE.parent.parent):
    sys.path.insert(0, str(candidate))
from qc import summarize  # noqa: E402

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()


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


def require(condition, message):
    if not condition:
        raise ValueError(message)


def available_kib():
    fields = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    return int(fields['MemAvailable'].split()[0])


def validate(plan):
    require(plan['format'] == 'brunello-full-run-v1', 'plan format drift')
    require(plan['authorization'].get('approved_by') and plan['authorization'].get('decision_ref'),
            'full-run plan requires an explicit recorded user authorization')
    require(plan['guide_count'] == 77441, 'full library only')
    require(digest(plan['guides']) == plan['pins'][plan['guides']], 'library checksum drift')
    m, b, z = plan['config']
    require([m, b, z] in [[2, 1, 1], [3, 1, 1], [4, 1, 1], [5, 1, 2]], 'unknown configuration')
    require(1 <= plan['chunk_size'] <= 77441, 'chunk size out of range')
    for assembly, spec in plan['assemblies'].items():
        for key in ['reference', 'binding', 'fai', 'scoring_reference_index']:
            require(Path(spec[key]).is_absolute() and spec[key] in plan['pins'], 'unpinned '+key)
        require(Path(spec['index']).is_absolute(), 'index prefix must be absolute')
        for suffix in SUFFIXES:
            require(spec['index']+'.'+suffix in plan['pins'], 'unpinned component '+suffix)
        binding = json.loads(Path(spec['binding']).read_text())
        require(binding['reference_sha256'] == spec['reference_sha256'], 'binding mismatch: '+assembly)
    return plan


def verify_pins(plan):
    for path, expected in plan['pins'].items():
        require(digest(path) == expected, 'SHA256 mismatch: '+path)


def run(command, out, err, metrics, plan, scratch, phase):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    with Path(out).open('x') as o, Path(err).open('x') as e:
        proc = subprocess.Popen(list(map(str, command)), stdout=o, stderr=e, start_new_session=True)
        try:
            while proc.poll() is None:
                require(time.monotonic()-start < plan['limits']['chunk_timeout_seconds'], phase+' timeout')
                require(available_kib() >= plan['limits']['min_running_available_kib'], phase+' available-memory guard')
                require(shutil.disk_usage(scratch).free >= plan['limits']['min_running_scratch_bytes']
                        and shutil.disk_usage(plan['root']).free >= plan['limits']['min_running_root_bytes'],
                        phase+' free-space guard')
                time.sleep(5)
            require(proc.returncode == 0, phase+' returned '+str(proc.returncode))
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    save(metrics, dict(phase=phase, wall_seconds=time.monotonic()-start, command=list(map(str, command)),
                      user_seconds=after.ru_utime-before.ru_utime, system_seconds=after.ru_stime-before.ru_stime,
                      cumulative_max_individual_descendant_rss_kib=after.ru_maxrss,
                      caveat='ru_maxrss is cumulative max individual child; tree peak is supervised.'))


def chunk_signature(plan, assembly, chunk, guides_sha):
    essential = {p: s for p, s in plan['pins'].items()
                 if p in [plan['cli'], plan['native']] or p.startswith(plan['assemblies'][assembly]['index']+'.')
                 or p == plan['assemblies'][assembly]['binding']}
    spec = dict(assembly=assembly, parameters=plan['config'], chunk=chunk,
                fraction=.75, pam='GG', chunk_guide_sha256=guides_sha, pins=essential)
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def inner(plan, assembly, attempt):
    verify_pins(plan)
    spec = plan['assemblies'][assembly]
    scratch_root = Path('/scratch')/os.environ['USER']
    scratch_root.mkdir(mode=0o700, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix='brunello-full-'+os.environ['SLURM_JOB_ID']+'-', dir=str(scratch_root)))
    m, b, z = plan['config']
    resource.setrlimit(resource.RLIMIT_FSIZE, (plan['limits']['file_limit_bytes'],)*2)
    os.environ.update(TMPDIR=str(scratch), PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',
                      OPENBLAS_NUM_THREADS='1', TF_NUM_INTRAOP_THREADS='1', TF_NUM_INTEROP_THREADS='1',
                      CUDA_VISIBLE_DEVICES='-1')
    status = attempt/'status.json'
    try:
        names = [line[1:] for line in Path(plan['guides']).read_text().splitlines() if line.startswith('>')]
        require(len(names) == plan['guide_count'], 'library record count drift')
        chunks_dir = attempt/'chunks'
        chunks_dir.mkdir(exist_ok=True)
        merged = attempt/'eligible.cfd.paf'
        with merged.open('x') as out:
            for lo in range(0, plan['guide_count'], plan['chunk_size']):
                hi = min(lo+plan['chunk_size'], plan['guide_count'])
                chunk = lo//plan['chunk_size']
                chunk_paf = chunks_dir/('chunk-%04d.paf' % chunk)
                checkpoint = chunks_dir/('chunk-%04d.complete.json' % chunk)
                chunk_guides = chunks_dir/('chunk-%04d.fa' % chunk)
                if not chunk_guides.is_file():
                    lines = Path(plan['guides']).read_text().splitlines()
                    with chunk_guides.open('x') as f:
                        for i in range(2*lo, 2*hi):
                            f.write(lines[i]+'\n')
                guides_sha = digest(chunk_guides)
                if checkpoint.is_file():
                    data = json.loads(checkpoint.read_text())
                    require(data['search_signature'] == chunk_signature(plan, assembly, chunk, guides_sha),
                            'checkpoint scientific identity drift')
                    require(digest(chunk_paf) == data['paf_sha256'], 'chunk output changed: '+str(chunk_paf))
                    with chunk_paf.open() as src, out:
                        shutil.copyfileobj(src, out)
                    continue
                command = [plan['cli'], '--search-backend', 'pam-anchored', '--pam-search-bin', plan['native'],
                           '--columba-index', spec['index'], '-r', spec['reference'], '-G', chunk_guides,
                           '-m', m, '-b', b, '-z', z, '-f', '.75', '-p', 'GG', '-t', 1,
                           '--mismatch-scores', plan['matrices']['mismatch'], '--pam-scores', plan['matrices']['pam']]
                run(command, chunk_paf, chunks_dir/('chunk-%04d.stderr' % chunk),
                    chunks_dir/('chunk-%04d.metrics.json' % chunk), plan, scratch, 'search-chunk-%d' % chunk)
                save(checkpoint, dict(search_signature=chunk_signature(plan, assembly, chunk, guides_sha),
                     paf=str(chunk_paf), paf_sha256=digest(chunk_paf),
                     rows=sum(1 for _ in chunk_paf.open())))
                save(status, dict(stage='search', completed_chunks=chunk+1,
                     total_chunks=(plan['guide_count']+plan['chunk_size']-1)//plan['chunk_size']))
                with chunk_paf.open() as src, out:
                    shutil.copyfileobj(src, out)
        audit_ref = scratch/'audit-reference.fa'
        audit_ref.symlink_to(spec['reference'])
        shutil.copyfile(spec['fai'], str(audit_ref)+'.fai')
        audit = attempt/'independent-audit'
        command = ['/usr/bin/python3', Path(plan['release_scripts'])/'audit_emitted.py', '--paf', merged,
                   '--fasta', audit_ref, '--guides', plan['guides'], '--matrices', Path(plan['matrices']['mismatch']).parent,
                   '--samtools', plan['samtools'], '--output', audit,
                   '--max-mismatches', m, '--max-bulges', b, '--max-bulge-size', z, '--min-match-fraction', '.75']
        run(command, attempt/'audit.stdout', attempt/'audit.stderr', attempt/'audit-metrics.json', plan, scratch, 'audit')
        audited = json.loads((audit/'summary.json').read_text())
        require(audited['status'] == 'PASS' and audited['errors'] == 0, 'independent audit failed')
        annotated = attempt/'eligible.both.paf'
        command = [plan['cli'], '-r', spec['reference'], '-G', plan['guides'], '--score-paf', merged,
                   '--scoring', 'both', '--scoring-output-dir', attempt/'neural',
                   '--crispr-bulge-python', plan['model_python'], '--crispr-bulge-model-root', plan['models'],
                   '--scoring-script', plan['scoring_script'],
                   '--scoring-reference-sha256', spec['reference_sha256'],
                   '--scoring-reference-index', spec['scoring_reference_index'],
                   '--scoring-reference-index-sha256', plan['pins'][spec['scoring_reference_index']],
                   '--scoring-worker-timeout', plan['limits']['scoring_timeout_seconds']]
        run(command, annotated, attempt/'scoring.stderr', attempt/'scoring-metrics.json', plan, scratch, 'scoring')
        require(digest(annotated) == digest(attempt/'neural/scoring/output.paf'), 'scoring bundle mismatch')
        result = summarize(merged, annotated, attempt/'neural/scoring/annotations.tsv', audit/'records.tsv', names)
        result.update(assembly=assembly, config=plan['config'], guide_count=plan['guide_count'])
        save(attempt/'qc-summary.json', result)
        verify_pins(plan)
        outputs = {str(p): dict(sha256=digest(p), bytes=p.stat().st_size)
                   for p in attempt.rglob('*') if p.is_file() and p.parent != chunks_dir}
        save(attempt/'STAGES_COMPLETE.json', dict(status='STAGES_COMPLETE_AWAITING_SUPERVISOR', outputs=outputs))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        for p in (attempt/'chunks').glob('chunk-*.paf'):
            p.unlink()  # merged PAF already contains every completed chunk; keep checkpoints/metrics only


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', required=True, type=Path)
    p.add_argument('--assembly', required=True)
    p.add_argument('--inner', type=Path)
    a = p.parse_args()
    plan = validate(json.loads(a.plan.read_text()))
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURMD_NODENAME'), 'Slurm execution required')
    if a.inner:
        inner(plan, a.assembly, a.inner)
        return
    logical = Path(plan['root'])/'02_search'/a.assembly
    logical.mkdir(parents=True, exist_ok=True)
    with (logical/'task.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        success = logical/'SUCCESS.json'
        if success.exists():
            result = json.loads(success.read_text())
            require(result['plan_sha256'] == digest(a.plan), 'completed task has a different plan')
            for path, spec in result['outputs'].items():
                require(digest(path) == spec['sha256'], 'completed output changed: '+path)
            print('Already complete; checksums verified; no rerun.')
            return
        attempt = logical/('attempt-'+os.environ['SLURM_JOB_ID'])
        attempt.mkdir()
        save(attempt/'INCOMPLETE.json', dict(plan=str(a.plan), assembly=a.assembly))
        try:
            verify_pins(plan)
            command = list(map(str, ['/usr/bin/python3', plan['supervisor'], '--limit-kib', plan['limits']['rss_kib'],
                       '--interval', '.1', '--log', attempt/'monitor.log', '--summary', attempt/'monitor-summary.json',
                       '--', '/usr/bin/python3', Path(__file__).resolve(), '--plan', a.plan,
                       '--assembly', a.assembly, '--inner', attempt]))
            subprocess.run(command, check=True)
            monitor = json.loads((attempt/'monitor-summary.json').read_text())
            require(monitor['reason'] == 'completed', 'supervisor did not pass')
            data = json.loads((attempt/'STAGES_COMPLETE.json').read_text())
            for path, spec in data['outputs'].items():
                require(digest(path) == spec['sha256'], 'output checksum mismatch')
            save(success, dict(status='COMPLETE', plan_sha256=digest(a.plan), attempt=str(attempt),
                               outputs=data['outputs'], monitor=monitor))
            (attempt/'INCOMPLETE.json').unlink()
        except BaseException as exc:
            save(attempt/'FAILED.json', dict(error=str(exc)))
            raise


if __name__ == '__main__':
    main()
