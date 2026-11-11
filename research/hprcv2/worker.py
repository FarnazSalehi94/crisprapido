#!/usr/bin/env python3
"""Exclusive, guarded PILOTS only. No entry point authorizes full Brunello."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time

from qc import summarize, require


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


def validate(plan):
    require(plan['format'] == 'hprcv2-guarded-pilot-v1', 'unsupported plan')
    require(plan['assembly'] == 'HG00735_pat', 'this safety pilot is restricted to the already bound assembly')
    require(plan['project'] in ('sbds', 'brunello'), 'unknown project')
    require(Path(plan['root']).name == plan['project']+'_pam_aware_v1', 'cross-project output root')
    require(plan['full_brunello_authorized'] is False, 'full Brunello is not authorized')
    require(plan['configurations'] == ({'primary': [5, 1, 2], 'historical-domain-control': [2, 1, 1]}
            if plan['project'] == 'sbds' else {'A': [4, 1, 1], 'B': [5, 1, 2]}), 'parameter drift')
    require(plan['fraction'] == .75 and plan['pam'] == 'GG', 'eligibility drift')
    require(plan['guide_count'] == (1 if plan['project'] == 'sbds' else 3), 'pilot guide count drift')
    require(plan['reference_sha256'] == 'bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b', 'reference identity drift')
    require(plan['limits'] == dict(rss_kib=50331648, min_start_available_kib=67108864,
            min_running_available_kib=8388608, min_scratch_bytes=60000000000,
            min_running_disk_bytes=5000000000, stage_timeout_seconds=6600, file_limit_bytes=10737418240), 'resource guard drift')
    for path, expected in plan['pins'].items():
        require(Path(path).is_absolute() and len(expected) == 64, 'missing absolute path/pin')
    require(str(Path(__file__).resolve()) in plan['pins'], 'worker source must be pinned')
    return plan


def verify_pins(plan):
    for path, expected in plan['pins'].items():
        require(digest(path) == expected, 'SHA256 mismatch: '+path)
    binding = json.loads(Path(plan['index']+'.crispr-binding.json').read_text())
    require(binding['reference_sha256'] == plan['reference_sha256'], 'native binding/reference mismatch')


def run(command, output, error, metrics, plan, scratch):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    with Path(output).open('x') as out, Path(error).open('x') as err:
        proc = subprocess.Popen(list(map(str, command)), stdout=out, stderr=err, start_new_session=True)
        try:
            while proc.poll() is None:
                require(time.monotonic()-start < plan['limits']['stage_timeout_seconds'], 'stage timeout')
                require(available_kib() >= plan['limits']['min_running_available_kib'], 'shared-node available-memory guard')
                require(min(shutil.disk_usage(scratch).free, shutil.disk_usage(plan['root']).free)
                        >= plan['limits']['min_running_disk_bytes'], 'free-space guard')
                time.sleep(1)
            require(proc.returncode == 0, 'stage returned '+str(proc.returncode)+': '+str(command[0]))
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    save(metrics, dict(command=list(map(str, command)), wall_seconds=time.monotonic()-start,
                      user_seconds=after.ru_utime-before.ru_utime, system_seconds=after.ru_stime-before.ru_stime,
                      cumulative_max_individual_descendant_rss_kib=after.ru_maxrss,
                      caveat='ru_maxrss is cumulative maximum individual child, not per-stage aggregate or a hard cap. Tree peak is measured by supervisor.'))


def search_signature(plan, task):
    essential = {path: sha for path, sha in plan['pins'].items()
                 if path in [plan['cli'], plan['native'], plan['index']+'.crispr-binding.json']
                 or path.startswith(plan['source']+'/')}
    spec = dict(assembly=plan['assembly'], project=plan['project'], parameters=plan['configurations'][task],
                fraction=plan['fraction'], pam=plan['pam'], guide_sha256=plan['pins'][plan['guides']],
                reference_sha256=plan['reference_sha256'], pins=essential)
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def inner(plan, task, attempt, resume_search=None):
    verify_pins(plan)
    limits = plan['limits']
    require(available_kib() >= limits['min_start_available_kib'], 'insufficient shared-node available memory at start')
    scratch_root = Path('/scratch')/os.environ['USER']
    scratch_root.mkdir(mode=0o700, exist_ok=True)
    require(shutil.disk_usage(scratch_root).free >= limits['min_scratch_bytes'], 'insufficient scratch')
    scratch = Path(tempfile.mkdtemp(prefix='hprcv2-'+os.environ['SLURM_JOB_ID']+'-', dir=str(scratch_root)))
    save(attempt/'preflight.json', dict(host=socket.gethostname(), job=os.environ['SLURM_JOB_ID'],
              scratch=str(scratch), mem_available_kib=available_kib(), scratch_free_bytes=shutil.disk_usage(scratch).free,
              plan=plan, native_search='PAM-aware Columba traversal only; no WFA2 candidate verifier; no CFD pruning.'))
    resource.setrlimit(resource.RLIMIT_FSIZE, (limits['file_limit_bytes'], limits['file_limit_bytes']))
    os.environ.update(TMPDIR=str(scratch), PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',
                      OPENBLAS_NUM_THREADS='1', TF_NUM_INTRAOP_THREADS='1', TF_NUM_INTEROP_THREADS='1',
                      CUDA_VISIBLE_DEVICES='-1')
    try:
        guides = attempt/'guides.fa'
        shutil.copyfile(plan['guides'], guides)
        require(digest(guides) == plan['pins'][plan['guides']], 'guide staging mismatch')
        lines = guides.read_text().splitlines()
        require(len(lines) == 2*plan['guide_count'], 'pilot must have exact two-line FASTA records')
        require(all(n.startswith('>') and len(s) == 20 and set(s) <= set('ACGT') for n, s in zip(lines[::2], lines[1::2])), 'invalid guide record')
        names = [line[1:] for line in lines[::2]]
        require(len(set(names)) == plan['guide_count'], 'guide count/ID mismatch before search')
        if plan['project'] == 'brunello':
            full = Path(plan['full_library']).read_text().splitlines()
            require(len(full) == 2*77441 and set(zip(lines[::2], lines[1::2])) <= set(zip(full[::2], full[1::2])), 'pilot not in frozen Brunello library')
        if plan['project'] == 'sbds':
            require(guides.read_text() == '>SBDS\nACTGAAATCTGTAAGCAGGC\n', 'SBDS input drift')
        score_root, qc_root = attempt, attempt
        if plan['project'] == 'sbds':
            score_root = Path(plan['root'])/'02_scoring'/plan['assembly']/task/attempt.name
            qc_root = Path(plan['root'])/'03_qc'/plan['assembly']/task/attempt.name
            score_root.mkdir(parents=True)
            qc_root.mkdir(parents=True)
        source = Path(plan['source'])
        audit_ref = scratch/'audit-reference.fa'
        audit_ref.symlink_to(plan['reference'])
        shutil.copyfile(plan['fai'], str(audit_ref)+'.fai')
        m, b, z = plan['configurations'][task]
        run([plan['model_python'], '-c', 'import json, tensorflow as tf, numpy as np, h5py; print(json.dumps(dict(tensorflow=tf.__version__, numpy=np.__version__, h5py=h5py.__version__)))'],
            attempt/'inference-environment.json', attempt/'environment.stderr', attempt/'environment-metrics.json', plan, scratch)
        command = [plan['cli'], '--search-backend', 'pam-anchored', '--pam-search-bin', plan['native'],
                   '--columba-index', plan['index'], '-r', plan['reference'], '-G', guides,
                   '-m', m, '-b', b, '-z', z, '-f', '.75', '-p', 'GG', '-t', 1,
                   '--mismatch-scores', source/'mismatch_scores.txt', '--pam-scores', source/'pam_scores.txt']
        raw = attempt/'eligible.cfd.paf'
        signature = search_signature(plan, task)
        if resume_search:
            checkpoint = json.loads(resume_search.read_text())
            require(checkpoint['status'] == 'SEARCH_PROCESS_COMPLETE_AWAITING_QC' and checkpoint['search_signature'] == signature,
                    'resume checkpoint has a different scientific search identity')
            previous = Path(checkpoint['paf'])
            shutil.copyfile(previous, raw)
            require(digest(raw) == checkpoint['paf_sha256'], 'resume PAF checksum mismatch')
            save(attempt/'search-metrics.json', dict(reused_search=True, original_checkpoint=str(resume_search),
                 original_checkpoint_sha256=digest(resume_search), original_search_metrics=checkpoint['metrics']))
        else:
            run(command, raw, attempt/'search.stderr', attempt/'search-metrics.json', plan, scratch)
        save(attempt/'SEARCH_COMPLETE.json', dict(status='SEARCH_PROCESS_COMPLETE_AWAITING_QC',
             search_signature=signature, paf=str(raw), paf_sha256=digest(raw),
             metrics=json.loads((attempt/'search-metrics.json').read_text())))
        audit_command = ['/usr/bin/python3', source/'scripts/release_validation/audit_emitted.py',
                '--paf', raw, '--fasta', audit_ref, '--guides', guides, '--matrices', source,
                '--samtools', plan['samtools'], '--output', qc_root/'independent-audit',
                '--max-mismatches', m, '--max-bulges', b, '--max-bulge-size', z, '--min-match-fraction', '.75']
        run(audit_command, qc_root/'audit.stdout', qc_root/'audit.stderr', qc_root/'audit-metrics.json', plan, scratch)
        command = [plan['cli'], '-r', plan['reference'], '-G', guides, '--score-paf', raw,
                '--scoring', 'both', '--scoring-output-dir', score_root/'neural',
                '--crispr-bulge-python', plan['model_python'], '--crispr-bulge-model-root', plan['models'],
                '--scoring-script', source/'research/crispr_bulge/score_paf.py',
                '--scoring-reference-sha256', plan['reference_sha256'],
                '--scoring-reference-index', plan['reference_index'],
                '--scoring-reference-index-sha256', plan['pins'][plan['reference_index']],
                '--scoring-worker-timeout', 1800]
        annotated = score_root/'eligible.both.paf'
        run(command, annotated, score_root/'scoring.stderr', score_root/'scoring-metrics.json', plan, scratch)
        require(digest(annotated) == digest(score_root/'neural/scoring/output.paf'), 'scoring bundle/stdout mismatch')
        result = summarize(raw, annotated, score_root/'neural/scoring/annotations.tsv', qc_root/'independent-audit/records.tsv', names)
        result.update(project=plan['project'], task=task, assembly=plan['assembly'], parameters=dict(m=m, b=b, z=z, f=.75, p='GG'),
                      plan_sha256=plan['_plan_sha256'])
        save(qc_root/'qc-summary.json', result)
        verify_pins(plan)
        outputs = {str(p): dict(sha256=digest(p), bytes=p.stat().st_size)
                   for directory in set([attempt, score_root, qc_root]) for p in directory.rglob('*')
                   if p.is_file() and p.name not in ['monitor.log', 'monitor-summary.json', 'INCOMPLETE.json']}
        save(attempt/'STAGES_COMPLETE.json', dict(status='STAGES_COMPLETE_AWAITING_SUPERVISOR', outputs=outputs,
                  total_stored_bytes=sum(r['bytes'] for r in outputs.values()), caveat='Includes duplicate diagnostic/scoring PAF and TSV; not a minimal-retention estimate.'))
    finally:
        shutil.rmtree(scratch)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', required=True, type=Path)
    p.add_argument('--task', required=True)
    p.add_argument('--inner', type=Path)
    p.add_argument('--resume-search', type=Path, help='Explicit checksum-bound completed search checkpoint; never guess from a partial PAF')
    a = p.parse_args()
    plan = validate(json.loads(a.plan.read_text()))
    plan['_plan_sha256'] = digest(a.plan)
    require(a.task in plan['configurations'], 'unknown task')
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURMD_NODENAME'), 'Slurm execution required')
    require(int(os.environ['SLURM_CPUS_PER_TASK']) >= 8, 'reserve at least eight CPUs; actual application uses one')
    if a.inner:
        inner(plan, a.task, a.inner, a.resume_search)
        return
    stage = '01_search' if plan['project'] == 'sbds' else '01_pilot'
    logical = Path(plan['root'])/stage/plan['assembly']/a.task
    logical.mkdir(parents=True, exist_ok=True)
    with (logical/'task.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        success = logical/'SUCCESS.json'
        if success.exists():
            result = json.loads(success.read_text())
            require(result['plan_sha256'] == plan['_plan_sha256'], 'completed task has a different plan; use new namespace')
            for path, spec in result['outputs'].items():
                require(digest(path) == spec['sha256'], 'completed output changed: '+path)
            print('Already complete; checksums verified; no experiment rerun.')
            return
        attempt = logical/('attempt-'+os.environ['SLURM_JOB_ID'])
        attempt.mkdir()
        save(attempt/'INCOMPLETE.json', dict(plan=str(a.plan), task=a.task))
        try:
            verify_pins(plan)
            command = ['/usr/bin/python3', plan['supervisor'], '--limit-kib', str(plan['limits']['rss_kib']),
                       '--interval', '.1', '--log', str(attempt/'monitor.log'), '--summary', str(attempt/'monitor-summary.json'),
                       '--', '/usr/bin/python3', str(Path(__file__).resolve()), '--plan', str(a.plan), '--task', a.task, '--inner', str(attempt)]
            if a.resume_search:
                command.extend(['--resume-search', str(a.resume_search.resolve())])
            subprocess.run(command, check=True)
            monitor = json.loads((attempt/'monitor-summary.json').read_text())
            require(monitor['reason'] == 'completed', 'monitor did not pass')
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
