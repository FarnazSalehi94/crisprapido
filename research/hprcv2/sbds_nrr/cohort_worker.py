#!/usr/bin/env python3
"""Restartable SBDS NGG/NRR per-assembly worker. No full-cohort entry point.

Per assembly, runs the four literal-PAM searches (m5/b1/z2/f0.75), an
independent audit with the explicit allowed PAM, and both-score annotation.
Searches are restart-safe: each PAM's final PAF lives at a stable path with a
checksum-verified SEARCH_COMPLETE checkpoint; a re-submitted job reuses verified
searches and re-runs only cheap audit/scoring stages. Uncommitted partial PAFs
from a crashed attempt are removed, never silently trusted. One CPU; the whole
job runs under the supervisor's memory/disk guards. Attempts are exclusive;
SUCCESS markers are checksum-verified on re-submission.
"""
import argparse
import collections
import csv
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import socket
import subprocess
import sys
import tempfile
import time

import signal

PAMS = ('NGG', 'NAG', 'NGA', 'NAA')
SUFFIX = dict(zip(PAMS, ('GG', 'AG', 'GA', 'AA')))


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def save(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2); f.write('\n')


def available_kib():
    fields = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    return int(fields['MemAvailable'].split()[0])


def require(cond, msg):
    if not cond: raise ValueError(msg)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod

def verify_pins(plan, extra=None):
    for path, expected in list(plan['pins'].items())+list((extra or {}).items()):
        require(digest(path) == expected, 'SHA256 mismatch: '+path)


def run(command, output, error, metrics, limits, guard_paths):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    with Path(output).open('x') as out, Path(error).open('x') as err:
        proc = subprocess.Popen(list(map(str, command)), stdout=out, stderr=err, start_new_session=True)
        try:
            while proc.poll() is None:
                require(time.monotonic()-start < limits['stage_timeout_seconds'], 'stage timeout')
                require(available_kib() >= limits['min_running_available_kib'], 'shared-node available-memory guard')
                for path in guard_paths:
                    require(shutil.disk_usage(str(path)).free >= limits['min_running_disk_bytes'],
                            'free-space guard: '+str(path))
                time.sleep(1)
            require(proc.returncode == 0, 'stage returned '+str(proc.returncode)+': '+str(command[0]))
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try: proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL); proc.wait()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    save(metrics, dict(command=list(map(str, command)), wall_seconds=time.monotonic()-start,
        user_seconds=after.ru_utime-before.ru_utime, system_seconds=after.ru_stime-before.ru_stime,
        cumulative_max_individual_descendant_rss_kib=after.ru_maxrss,
        caveat='ru_maxrss is cumulative maximum individual child; supervisor measures tree peak.'))


def search_signature(plan, entry, suffix):
    essential = {path: sha for path, sha in plan['pins'].items()
                 if path in [entry['cli'], entry['native'], entry['index']+'.crispr-binding.json']
                 or path.startswith(plan['source']+'/')}
    return hashlib.sha256(json.dumps(dict(assembly=entry['assembly'],
        guide_sha256=plan['pins'][plan['guides']], parameters=[5, 1, 2, .75], pam=suffix,
        reference_sha256=entry['reference_sha256'], pins=essential), sort_keys=True).encode()).hexdigest()


def inner(plan, assembly, attempt, scratch):
    entry = plan['assemblies'][assembly]
    limits = plan['limits']
    root = Path(plan['root'])
    verify_pins(plan, entry['pins'])
    os.environ.update(TMPDIR=str(scratch), PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',
                      OPENBLAS_NUM_THREADS='1', TF_NUM_INTRAOP_THREADS='1', TF_NUM_INTEROP_THREADS='1',
                      CUDA_VISIBLE_DEVICES='-1')
    resource.setrlimit(resource.RLIMIT_FSIZE, (limits['file_limit_bytes'],)*2)
    source = Path(plan['source'])
    save(attempt/'preflight.json', dict(host=socket.gethostname(), job=os.environ['SLURM_JOB_ID'],
        scratch=str(scratch), mem_available_kib=available_kib(), plan_sha256=plan['_plan_sha256']))
    guards = [scratch, root]
    auditor = load('cohort_auditor', plan['auditor'])
    for pam in PAMS:
        suffix = SUFFIX[pam]
        logical = root/'02_search'/assembly/pam
        logical.mkdir(parents=True, exist_ok=True)
        paf = logical/'eligible.cfd.paf'
        checkpoint = logical/'SEARCH_COMPLETE.json'
        signature = search_signature(plan, entry, suffix)
        if checkpoint.exists():
            data = json.loads(checkpoint.read_text())
            require(data['search_signature'] == signature, 'search checkpoint signature drift for '+pam)
            require(digest(paf) == data['paf_sha256'] == digest(data['paf']) if paf.exists() else False,
                    'search checkpoint PAF mismatch for '+pam)
            save(attempt/(pam+'.search-metrics.json'),
                 dict(reused_search=True, checkpoint=str(checkpoint),
                      original_search_metrics=data.get('metrics')))
        else:
            if paf.exists():
                paf.unlink()  # uncommitted partial from a crashed attempt
            run([entry['cli'], '--search-backend', 'pam-anchored', '--pam-search-bin', entry['native'],
                 '--columba-index', entry['index'], '-r', entry['reference'], '-G', plan['guides'],
                 '-m', 5, '-b', 1, '-z', 2, '-f', '.75', '-p', suffix, '-t', 1,
                 '--mismatch-scores', source/'mismatch_scores.txt',
                 '--pam-scores', source/'pam_scores.txt'],
                paf, attempt/(pam+'.search.stderr'), attempt/(pam+'.search-metrics.json'), limits, guards)
            save(checkpoint, dict(search_signature=signature, paf=str(paf), paf_sha256=digest(paf),
                 metrics=json.loads((attempt/(pam+'.search-metrics.json')).read_text())))
        ref = scratch/('audit-'+pam+'.fa')
        if ref.exists() or ref.is_symlink(): ref.unlink()
        ref.symlink_to(entry['reference'])
        shutil.copyfile(entry['fai'], str(ref)+'.fai')
        audit_out = root/'04_qc'/assembly/pam/attempt.name
        if audit_out.exists(): shutil.rmtree(audit_out)
        audit_out.parent.mkdir(parents=True, exist_ok=True)
        audit = auditor.audit(paf, ref, plan['guides'], source, plan['samtools'], audit_out,
                              max_mismatches=5, max_bulges=1, max_bulge_size=2, allowed_pams=(suffix,))
        require(audit['status'] == 'PASS', 'independent audit failed for '+pam)
        score_out = root/'03_scoring'/assembly/pam/attempt.name
        if score_out.exists(): shutil.rmtree(score_out)
        score_out.mkdir(parents=True)
        run([entry['cli'], '-r', entry['reference'], '-G', plan['guides'], '--score-paf', paf,
             '--scoring', 'both', '--scoring-output-dir', score_out/'neural',
             '--crispr-bulge-python', plan['model_python'], '--crispr-bulge-model-root', plan['models'],
             '--scoring-script', plan.get('scoring_script') or source/'research/crispr_bulge/score_paf.py',
             '--scoring-reference-sha256', entry['reference_sha256'],
             '--scoring-reference-index', entry['reference_index'],
             '--scoring-reference-index-sha256', entry['pins'][entry['reference_index']],
             '--scoring-worker-timeout', 1800],
            score_out/'eligible.both.paf', score_out/'scoring.stderr', score_out/'scoring-metrics.json', limits, guards)
        require(digest(score_out/'eligible.both.paf') == digest(score_out/'neural/scoring/output.paf'),
                'scoring bundle/stdout mismatch')
        counts = collections.Counter()
        with (score_out/'neural/scoring/annotations.tsv').open() as f:
            for r in csv.DictReader(f, delimiter='\t'):
                require(r['pam'][1:] == suffix, 'annotation PAM drift')
                counts['rows'] += 1
                counts['cfd_defined'] += r['cfd'] != 'N/A'
                counts['cb_'+r['crispr_bulge_reason']] += 1
        require(counts['rows'] == audit['rows'], 'scoring row loss')
        save(attempt/(pam+'.summary.json'), dict(counts=dict(counts), audit_status=audit['status']))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--assembly', required=True)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    plan['_plan_sha256'] = digest(a.plan)
    require(plan['format'] == 'sbds-nrr-cohort-v1', 'unsupported plan')
    require(Path(plan['root']).name == 'sbds_ngg_nrr_v1', 'cross-project root')
    require(a.assembly in plan['assemblies'], 'assembly not in frozen plan')
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_CPUS_PER_TASK') == '1', 'one-CPU Slurm required')
    entry = plan['assemblies'][a.assembly]
    require(Path(str(entry['index'])+'.crispr-binding.json').exists(), 'binding sidecar missing')
    logical = Path(plan['root'])/'02_search'/a.assembly
    logical.mkdir(parents=True, exist_ok=True)
    with (logical/'assembly.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        success = logical/'SUCCESS.json'
        if success.exists():
            data = json.loads(success.read_text())
            for path, spec in data['outputs'].items():
                require(digest(path) == spec['sha256'], 'completed output changed: '+path)
            # Incremental expansion: a completed assembly is skipped under a NEW
            # plan only when its scientific search identity is unchanged.
            stored = data.get('search_signatures')
            if data['plan_sha256'] == plan['_plan_sha256']:
                pass
            elif stored:
                current = {pam: search_signature(plan, entry, SUFFIX[pam]) for pam in PAMS}
                require(stored == current, 'completed assembly has a different scientific search identity')
            else:
                raise ValueError('completed assembly was produced by a different plan and carries no '
                                 'recorded search signatures; refusing silent reuse')
            print('Already complete; checksums verified; no rerun.')
            return
        attempt = logical/('attempt-'+os.environ['SLURM_JOB_ID'])
        require(not attempt.exists(), 'attempt directory already exists for this job id')
        attempt.mkdir()
        save(attempt/'INCOMPLETE.json', dict(plan_sha256=plan['_plan_sha256']))
        scratch = None
        try:
            verify_pins(plan, plan['assemblies'][a.assembly]['pins'])
            limits = plan['limits']
            require(available_kib() >= limits['min_start_available_kib'], 'insufficient available memory at start')
            base = Path('/scratch')/os.environ['USER']
            base.mkdir(mode=0o700, exist_ok=True)
            require(shutil.disk_usage(base).free >= limits['min_scratch_bytes'], 'insufficient scratch')
            scratch = Path(tempfile.mkdtemp(prefix='sbds-nrr-'+os.environ['SLURM_JOB_ID']+'-', dir=base))
            inner(plan, a.assembly, attempt, scratch)
            verify_pins(plan, plan['assemblies'][a.assembly]['pins'])
            job = 'attempt-'+os.environ['SLURM_JOB_ID']
            outputs = {}
            for f in logical.rglob('*'):
                if (f.is_file() and f.name not in ('assembly.lock', 'INCOMPLETE.json', 'SUCCESS')
                        and 'attempt-' in str(f) and job not in str(f)):
                    continue
                if f.is_file() and f.name not in ('assembly.lock', 'INCOMPLETE.json', 'SUCCESS.json'):
                    outputs[str(f)] = dict(sha256=digest(f), bytes=f.stat().st_size)
            for stage in ('03_scoring', '04_qc'):
                for f in (Path(plan['root'])/stage/a.assembly).rglob('*'):
                    if f.is_file() and job in str(f):
                        outputs[str(f)] = dict(sha256=digest(f), bytes=f.stat().st_size)
            save(success, dict(status='COMPLETE', plan_sha256=plan['_plan_sha256'], attempt=str(attempt),
                 assembly=a.assembly, outputs=outputs, pams=PAMS,
                 search_signatures={pam: search_signature(plan, entry, SUFFIX[pam]) for pam in PAMS}))
            (attempt/'INCOMPLETE.json').unlink()
        except BaseException as exc:
            save(attempt/'FAILED.json', dict(error=str(exc), plan_sha256=plan['_plan_sha256']))
            raise
        finally:
            if scratch is not None:
                shutil.rmtree(scratch, ignore_errors=True)


if __name__ == '__main__':
    main()
