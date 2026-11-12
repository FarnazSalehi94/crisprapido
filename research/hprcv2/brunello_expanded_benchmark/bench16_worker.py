#!/usr/bin/env python3
"""Guarded 16-guide x 4-configuration x 2-assembly Brunello benchmark.

One Slurm task per (assembly, configuration): a single search invocation for
all guides, then independent audit, then optional separate neural annotation.
Restart-safe: a checksum-bound SEARCH_COMPLETE checkpoint allows explicit resume
of post-search stages without repeating a completed search. No full-library or
cohort launch is possible from this worker.
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

CONFIGS = dict(historical=[2, 1, 1], intermediate=[3, 1, 1],
               **{'candidate-A': [4, 1, 1], 'candidate-B': [5, 1, 2]})
LIMITS = dict(stage_timeout_seconds=10800, min_start_available_kib=67108864,
              min_running_available_kib=8388608, min_scratch_bytes=60000000000,
              min_running_disk_bytes=10000000000, file_limit_bytes=21474836480,
              rss_kib=50331648)
GUIDE_COUNT = 16
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
    require(plan['format'] == 'brunello-expanded-benchmark-v1', 'plan format drift')
    require(plan['limits'] == LIMITS, 'resource-guard drift')
    require(plan['configs'] == CONFIGS, 'configuration drift')
    require(plan['guide_count'] == GUIDE_COUNT, 'panel size drift')
    require(set(plan['assemblies']) == {'HG00735_pat', 'HG01496_pat'}, 'assembly set drift')
    require(plan['full_run_authorized'] is False, 'full-run authorization must stay false')
    for assembly, spec in plan['assemblies'].items():
        for key in ['reference', 'binding', 'fai', 'scoring_reference_index']:
            require(Path(spec[key]).is_absolute() and spec[key] in plan['pins'], 'unpinned '+key+' for '+assembly)
        # The index prefix itself is not a file; every component must be pinned instead.
        require(Path(spec['index']).is_absolute(), 'index prefix must be absolute: '+assembly)
        for suffix in SUFFIXES:
            require(spec['index']+'.'+suffix in plan['pins'], 'unpinned index component '+suffix+' for '+assembly)
        binding = json.loads(Path(spec['binding']).read_text())
        require(binding['reference_sha256'] == spec['reference_sha256'], 'binding/reference mismatch: '+assembly)
    require(len({a['reference_sha256'] for a in plan['assemblies'].values()}) == 2, 'assembly references must differ')
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
                require(time.monotonic()-start < plan['limits']['stage_timeout_seconds'], phase+' timeout')
                require(available_kib() >= plan['limits']['min_running_available_kib'], phase+' available-memory guard')
                require(min(shutil.disk_usage(scratch).free, shutil.disk_usage(plan['root']).free)
                        >= plan['limits']['min_running_disk_bytes'], phase+' free-space guard')
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
    save(metrics, dict(phase=phase, command=list(map(str, command)), wall_seconds=time.monotonic()-start,
                      user_seconds=after.ru_utime-before.ru_utime, system_seconds=after.ru_stime-before.ru_stime,
                      cumulative_max_individual_descendant_rss_kib=after.ru_maxrss,
                      caveat='ru_maxrss is cumulative maximum individual child, not aggregate; tree peak is supervised.'))


def parse_native(text):
    found = re.findall(r'^states=(\d+) search_wall_seconds=([0-9.eE+-]+)$', text, re.M)
    require(len(found) == 1, 'missing/duplicate native timing line')
    return int(found[0][0]), float(found[0][1])


def search_signature(plan, assembly, config):
    essential = {p: s for p, s in plan['pins'].items()
                 if p in [plan['cli'], plan['native']] or p.startswith(plan['assemblies'][assembly]['index']+'.')
                 or p == plan['assemblies'][assembly]['binding']}
    spec = dict(assembly=assembly, parameters=CONFIGS[config], fraction=.75, pam='GG',
                guide_sha256=plan['pins'][plan['guides']], pins=essential)
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def inner(plan, assembly, config, attempt, resume=None):
    verify_pins(plan)
    spec = plan['assemblies'][assembly]
    limits = plan['limits']
    require(available_kib() >= limits['min_start_available_kib'], 'insufficient available memory at start')
    scratch_root = Path('/scratch')/os.environ['USER']
    scratch_root.mkdir(mode=0o700, exist_ok=True)
    require(shutil.disk_usage(scratch_root).free >= limits['min_scratch_bytes'], 'insufficient scratch')
    scratch = Path(tempfile.mkdtemp(prefix='brunello16-'+os.environ['SLURM_JOB_ID']+'-', dir=str(scratch_root)))
    m, b, z = CONFIGS[config]
    save(attempt/'preflight.json', dict(host=socket.gethostname(), job=os.environ['SLURM_JOB_ID'],
          scratch=str(scratch), mem_available_kib=available_kib(), assembly=assembly, config=config,
          parameters=dict(m=m, b=b, z=z, f=.75, p='GG')))
    resource.setrlimit(resource.RLIMIT_FSIZE, (limits['file_limit_bytes'], limits['file_limit_bytes']))
    os.environ.update(TMPDIR=str(scratch), PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',
                      OPENBLAS_NUM_THREADS='1', TF_NUM_INTRAOP_THREADS='1', TF_NUM_INTEROP_THREADS='1',
                      CUDA_VISIBLE_DEVICES='-1')
    try:
        guides = attempt/'guides.fa'
        shutil.copyfile(plan['guides'], guides)
        require(digest(guides) == plan['pins'][plan['guides']], 'guide staging mismatch')
        audit_ref = scratch/'audit-reference.fa'
        audit_ref.symlink_to(spec['reference'])
        shutil.copyfile(spec['fai'], str(audit_ref)+'.fai')
        raw = attempt/'eligible.cfd.paf'
        signature = search_signature(plan, assembly, config)
        if resume is not None:
            checkpoint = json.loads(Path(resume).read_text())
            require(checkpoint['status'] == 'SEARCH_COMPLETE_AWAITING_QC' and checkpoint['search_signature'] == signature,
                    'resume checkpoint has a different scientific search identity')
            shutil.copyfile(checkpoint['paf'], raw)
            require(digest(raw) == checkpoint['paf_sha256'], 'resume PAF checksum mismatch')
            states, traversal = checkpoint['native_states'], checkpoint['native_traversal_seconds']
            save(attempt/'search-metrics.json', dict(reused_search=True, original_checkpoint=str(resume),
                 original_checkpoint_sha256=digest(Path(resume)), original_search_metrics=checkpoint['metrics']))
        else:
            command = [plan['cli'], '--search-backend', 'pam-anchored', '--pam-search-bin', plan['native'],
                       '--columba-index', spec['index'], '-r', spec['reference'], '-G', guides,
                       '-m', m, '-b', b, '-z', z, '-f', '.75', '-p', 'GG', '-t', 1,
                       '--mismatch-scores', plan['matrices']['mismatch'], '--pam-scores', plan['matrices']['pam']]
            run(command, raw, attempt/'search.stderr', attempt/'search-metrics.json', plan, scratch, 'search')
            states, traversal = parse_native((attempt/'search.stderr').read_text())
            metrics = json.loads((attempt/'search-metrics.json').read_text())
            save(attempt/'SEARCH_COMPLETE.json', dict(status='SEARCH_COMPLETE_AWAITING_QC',
                 search_signature=signature, paf=str(raw), paf_sha256=digest(raw), metrics=metrics,
                 native_states=states, native_traversal_seconds=traversal))
        audit = attempt/'independent-audit'
        command = ['/usr/bin/python3', Path(plan['release_scripts'])/'audit_emitted.py', '--paf', raw,
                   '--fasta', audit_ref, '--guides', guides, '--matrices', Path(plan['matrices']['mismatch']).parent,
                   '--samtools', plan['samtools'], '--output', audit,
                   '--max-mismatches', m, '--max-bulges', b, '--max-bulge-size', z, '--min-match-fraction', '.75']
        run(command, attempt/'audit.stdout', attempt/'audit.stderr', attempt/'audit-metrics.json', plan, scratch, 'audit')
        audited = json.loads((audit/'summary.json').read_text())
        require(audited['status'] == 'PASS' and audited['errors'] == 0, 'independent audit failed')
        score_root = attempt/'neural'
        command = [plan['cli'], '-r', spec['reference'], '-G', guides, '--score-paf', raw,
                   '--scoring', 'both', '--scoring-output-dir', score_root,
                   '--crispr-bulge-python', plan['model_python'], '--crispr-bulge-model-root', plan['models'],
                   '--scoring-script', plan['scoring_script'],
                   '--scoring-reference-sha256', spec['reference_sha256'],
                   '--scoring-reference-index', spec['scoring_reference_index'],
                   '--scoring-reference-index-sha256', plan['pins'][spec['scoring_reference_index']],
                   '--scoring-worker-timeout', 7200]
        annotated = attempt/'eligible.both.paf'
        run(command, annotated, attempt/'scoring.stderr', attempt/'scoring-metrics.json', plan, scratch, 'scoring')
        require(digest(annotated) == digest(score_root/'scoring/output.paf'), 'scoring bundle/stdout mismatch')
        names = [line[1:] for line in guides.read_text().splitlines() if line.startswith('>')]
        require(len(names) == plan['guide_count'], 'guide count mismatch')
        result = summarize(raw, annotated, score_root/'scoring/annotations.tsv', audit/'records.tsv', names)
        result.update(assembly=assembly, config=config, parameters=dict(m=m, b=b, z=z, f=.75, p='GG'),
                      native_states=states, native_traversal_seconds=traversal,
                      native_mean_seconds_per_guide=traversal/plan['guide_count'],
                      guide_count=plan['guide_count'], plan_sha256=plan['_plan_sha256'])
        save(attempt/'qc-summary.json', result)
        verify_pins(plan)
        outputs = {str(p): dict(sha256=digest(p), bytes=p.stat().st_size)
                   for p in attempt.rglob('*') if p.is_file()
                   and p.name not in ['monitor.log', 'monitor-summary.json', 'INCOMPLETE.json', 'FAILED.json']}
        save(attempt/'STAGES_COMPLETE.json', dict(status='STAGES_COMPLETE_AWAITING_SUPERVISOR', outputs=outputs))
    finally:
        shutil.rmtree(scratch)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', required=True, type=Path)
    p.add_argument('--assembly', required=True)
    p.add_argument('--config', required=True)
    p.add_argument('--inner', type=Path)
    p.add_argument('--resume-search', type=Path)
    a = p.parse_args()
    plan = validate(json.loads(a.plan.read_text()))
    plan['_plan_sha256'] = digest(a.plan)
    require(a.assembly in plan['assemblies'] and a.config in CONFIGS, 'unknown assembly/config')
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURMD_NODENAME'), 'Slurm execution required')
    require(int(os.environ['SLURM_CPUS_PER_TASK']) >= 8, 'reserve at least eight CPUs; application uses one')
    if a.inner:
        inner(plan, a.assembly, a.config, a.inner, a.resume_search)
        return
    logical = Path(plan['root'])/'01_pilot'/'expanded-benchmark-v1'/a.assembly/a.config
    logical.mkdir(parents=True, exist_ok=True)
    with (logical/'task.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        success = logical/'SUCCESS.json'
        if success.exists():
            result = json.loads(success.read_text())
            require(result['plan_sha256'] == plan['_plan_sha256'], 'completed task has a different plan')
            for path, spec in result['outputs'].items():
                require(digest(path) == spec['sha256'], 'completed output changed: '+path)
            print('Already complete; checksums verified; no benchmark rerun.')
            return
        attempt = logical/('attempt-'+os.environ['SLURM_JOB_ID'])
        attempt.mkdir()
        save(attempt/'INCOMPLETE.json', dict(plan=str(a.plan), assembly=a.assembly, config=a.config))
        try:
            verify_pins(plan)
            command = list(map(str, ['/usr/bin/python3', plan['supervisor'], '--limit-kib', plan['limits']['rss_kib'],
                       '--interval', '.1', '--log', attempt/'monitor.log', '--summary', attempt/'monitor-summary.json',
                       '--', '/usr/bin/python3', Path(__file__).resolve(), '--plan', a.plan,
                       '--assembly', a.assembly, '--config', a.config, '--inner', attempt]))
            if a.resume_search:
                command.extend(['--resume-search', str(a.resume_search.resolve())])
            subprocess.run(command, check=True)
            monitor = json.loads((attempt/'monitor-summary.json').read_text())
            require(monitor['reason'] == 'completed', 'supervisor did not pass')
            data = json.loads((attempt/'STAGES_COMPLETE.json').read_text())
            for path, spec in data['outputs'].items():
                require(digest(path) == spec['sha256'], 'output checksum mismatch')
            for name in ['monitor.log', 'monitor-summary.json', 'STAGES_COMPLETE.json']:
                path = attempt/name
                data['outputs'][str(path)] = dict(sha256=digest(path), bytes=path.stat().st_size)
            save(success, dict(status='COMPLETE', plan_sha256=plan['_plan_sha256'], attempt=str(attempt),
                               outputs=data['outputs'], monitor=monitor))
            (attempt/'INCOMPLETE.json').unlink()
        except BaseException as exc:
            save(attempt/'FAILED.json', dict(error=str(exc), plan_sha256=plan['_plan_sha256']))
            raise


if __name__ == '__main__':
    main()
