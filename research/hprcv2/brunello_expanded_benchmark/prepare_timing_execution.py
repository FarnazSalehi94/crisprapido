#!/usr/bin/env python3
"""Freeze the bounded timing worker and Slurm script; never submits a job."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def option(command, flag):
    if command.count(flag) != 1: raise ValueError('missing/repeated flag: '+flag)
    return command[command.index(flag)+1]


def prepare(source, dest, authorization):
    if not authorization: raise ValueError('explicit benchmark authorization required')
    original = json.loads(source.read_text())
    if original['guide_count'] != 10000 or original['full_cohort_authorized']:
        raise ValueError('not a bounded 10k plan')
    dest.mkdir(exist_ok=False)
    repo = Path(__file__).resolve().parents[3]
    supervisor = Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/00_manifest/nrr-tooling-v1/supervise.py')
    for src, name in [(Path(__file__).with_name('timing_worker.py'), 'timing_worker.py'),
                      (repo/'research/hprcv2/sbds_nrr/wave_guard.py', 'wave_guard.py'),
                      (supervisor, 'supervise.py')]:
        shutil.copyfile(src, dest/name)
    plan = copy.deepcopy(original)
    plan.update(format='brunello-timing-10000-v1', authorization=authorization,
                original_plan=str(source), original_plan_sha256=digest(source))
    ref = option(plan['commands']['search'], '-r')
    scorer = Path(option(plan['commands']['score'], '--scoring-script')).parent
    for p in [source, Path(ref+'.fai'), scorer/'model_manifest.json',
              Path(option(plan['commands']['score'], '--crispr-bulge-python')),
              dest/'timing_worker.py', dest/'wave_guard.py', dest/'supervise.py']:
        plan['pins'][str(p)] = digest(p)
    # The audit interpreter is a node-provided system python3; hashing it is
    # node-dependent, so it is recorded at runtime rather than pinned. The
    # audited script, samtools and every scientific input remain pinned.
    plan['recorded_not_pinned'] = plan['commands']['audit'][0]
    plan_path = dest/'execution-plan.json'
    plan_path.write_text(json.dumps(plan, indent=2)+'\n')
    output = Path(plan['output'])
    script = f'''#!/bin/bash
#SBATCH --job-name=brunello10k-timing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=72:00:00
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TF_NUM_INTRAOP_THREADS=1 TF_NUM_INTEROP_THREADS=1 CUDA_VISIBLE_DEVICES=-1
D={dest}
O={output}
(cd "$D" && sha256sum --check tooling.sha256)
# A failed attempt has its own guard record; preserve it on any later restart.
mkdir "$O/guard-$SLURM_JOB_ID"
python3 "$D/wave_guard.py" --wave-dir "$O/guard-$SLURM_JOB_ID" -- python3 "$D/supervise.py" --limit-kib 50331648 --interval .5 --log "$O/supervisor-$SLURM_JOB_ID.jsonl" --summary "$O/supervisor-$SLURM_JOB_ID.json" -- python3 "$D/timing_worker.py" --plan "$D/execution-plan.json"
python3 "$D/timing_worker.py" --plan "$D/execution-plan.json" --finalize "$O/supervisor-$SLURM_JOB_ID.json"
'''
    (dest/'submit.slurm').write_text(script)
    (dest/'tooling.sha256').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in sorted(dest.iterdir()) if p.is_file()))
    return dict(plan=str(plan_path), sha256=digest(plan_path), output=str(output), status='PREPARED_NOT_SUBMITTED')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--source', type=Path, required=True)
    p.add_argument('--destination', type=Path, required=True); p.add_argument('--authorization', required=True)
    a = p.parse_args(); print(json.dumps(prepare(a.source, a.destination, a.authorization)))
