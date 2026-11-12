#!/usr/bin/env python3
"""Bounded literal-PAM pilot; deliberately has no cohort-launch entry point."""
import argparse
import collections
import csv
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import tempfile


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--pam', choices=['GG', 'AG', 'GA', 'AA'], required=True)
    a = ap.parse_args()
    p = json.loads(a.plan.read_text())
    sys.path.insert(0, str(Path(p['worker']).parent))
    w = load('guarded_worker', p['worker'])
    w.verify_pins(p)
    w.require(p['assembly'] == 'HG00735_pat' and Path(p['root']).name == 'sbds_ngg_nrr_v1', 'scope drift')
    w.require(Path(p['guides']).read_text() == '>SBDS\nACTGAAATCTGTAAGCAGGC\n', 'guide drift')
    w.require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_CPUS_PER_TASK') == '1', 'one-CPU Slurm pilot required')
    logical = Path(p['root'])/'01_pilot'/('N'+a.pam)
    logical.mkdir(exist_ok=True)
    with (logical/'lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        success = logical/'SUCCESS.json'
        if success.exists():
            s = json.loads(success.read_text())
            w.require(s['plan_sha256'] == w.digest(a.plan), 'plan changed')
            for f, h in s['outputs'].items():
                w.require(w.digest(f) == h, 'completed output changed: '+f)
            return
        out = logical/('attempt-'+os.environ['SLURM_JOB_ID'])
        out.mkdir()
        scratch = None
        try:
            w.require(w.available_kib() >= p['limits']['min_start_available_kib'], 'available memory preflight')
            base = Path('/scratch')/os.environ['USER']
            base.mkdir(exist_ok=True)
            w.require(shutil.disk_usage(base).free >= p['limits']['min_scratch_bytes'], 'scratch preflight')
            scratch = Path(tempfile.mkdtemp(prefix='sbds-nrr-', dir=base))
            os.environ.update(TMPDIR=str(scratch), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', TF_NUM_INTRAOP_THREADS='1', TF_NUM_INTEROP_THREADS='1', CUDA_VISIBLE_DEVICES='-1')
            resource.setrlimit(resource.RLIMIT_FSIZE, (p['limits']['file_limit_bytes'],)*2)
            def run(cmd, stem):
                w.run(cmd, out/(stem+'.stdout'), out/(stem+'.stderr'), out/(stem+'.metrics.json'), p, scratch)
            source = Path(p['source'])
            run([p['cli'], '--search-backend', 'pam-anchored', '--pam-search-bin', p['native'], '--columba-index', p['index'], '-r', p['reference'], '-G', p['guides'], '-m', 5, '-b', 1, '-z', 2, '-f', '.75', '-p', a.pam, '-t', 1, '--mismatch-scores', source/'mismatch_scores.txt', '--pam-scores', source/'pam_scores.txt'], 'search')
            ref = scratch/'reference.fa'
            ref.symlink_to(p['reference'])
            shutil.copyfile(p['fai'], str(ref)+'.fai')
            auditor = load('nrr_auditor', p['auditor'])
            audit = auditor.audit(out/'search.stdout', ref, p['guides'], source, p['samtools'], out/'audit', max_mismatches=5, max_bulge_size=2, allowed_pams=(a.pam,))
            w.require(audit['status'] == 'PASS', 'independent audit failed')
            run([p['cli'], '-r', p['reference'], '-G', p['guides'], '--score-paf', out/'search.stdout', '--scoring', 'both', '--scoring-output-dir', out/'neural', '--crispr-bulge-python', p['model_python'], '--crispr-bulge-model-root', p['models'], '--scoring-script', source/'research/crispr_bulge/score_paf.py', '--scoring-reference-sha256', p['reference_sha256'], '--scoring-reference-index', p['reference_index'], '--scoring-reference-index-sha256', p['pins'][p['reference_index']], '--scoring-worker-timeout', 1800], 'scoring')
            counts = collections.Counter()
            keys = set()
            with (out/'neural/scoring/annotations.tsv').open() as f:
                for r in csv.DictReader(f, delimiter='\t'):
                    w.require(r['pam'][1:] == a.pam, 'annotation PAM drift')
                    key = (r['guide_id'], r['contig'], r['strand'], r['end_0based'] if r['strand'] == '+' else r['start_0based'])
                    w.require(key not in keys, 'duplicate anchor')
                    keys.add(key)
                    counts['rows'] += 1
                    counts['cfd_defined'] += r['cfd'] != 'N/A'
                    counts['cb_'+r['crispr_bulge_reason']] += 1
            w.require(counts['rows'] == audit['rows'], 'scoring row loss')
            w.verify_pins(p)
            w.save(out/'summary.json', dict(pam='N'+a.pam, counts=dict(counts), audit=audit['status'], neural_noncanonical_policy='N/A unsupported_pam; independent model acceptance/validity investigation pending'))
            outputs = {str(f): w.digest(f) for f in out.rglob('*') if f.is_file()}
            w.save(success, dict(status='PILOT_STAGES_COMPLETE_NOT_COHORT_APPROVAL', plan_sha256=w.digest(a.plan), outputs=outputs))
        except BaseException as exc:
            w.save(out/'FAILED.json', dict(error=str(exc)))
            raise
        finally:
            if scratch is not None:
                shutil.rmtree(scratch)


if __name__ == '__main__':
    main()
