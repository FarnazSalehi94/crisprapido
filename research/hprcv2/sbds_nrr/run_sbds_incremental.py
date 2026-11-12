#!/usr/bin/env python3
"""Incremental SBDS NGG/NRR campaign driver — PREPARES, never auto-launches.

Reads the authoritative cohort manifest, selects assemblies that are
READY_BOUND and not yet COMPLETE, freezes a fresh cohort plan for exactly
those assemblies (per-assembly search identity unchanged ⇒ completed ones are
skipped under any later expanded plan), writes the ASSEMBLIES_FILE, and
prints the exact sbatch command. Launching requires passing --submit AND a
verbatim authorization string; without both it is a dry run.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', type=Path, default=ROOT/'00_manifest/authoritative-manifest-v2.json')
    p.add_argument('--repo', type=Path, default=Path('/export/local/home/salehi/crisprapido-pam-validation'))
    p.add_argument('--plan-tag', required=True)
    p.add_argument('--concurrency', type=int, default=16)
    p.add_argument('--submit', action='store_true')
    p.add_argument('--authorization', default='')
    p.add_argument('--time-limit', default='12:00:00')
    a = p.parse_args()
    manifest = json.loads(a.manifest.read_text())
    assert manifest['cohort_size'] == 466 and manifest['mutually_exclusive']
    ready = [r['assembly'] for r in manifest['rows'] if r['status'] == 'READY_BOUND']
    pending = [x for x in ready if not (ROOT/'02_search'/x/'SUCCESS.json').exists()]
    if not pending:
        print(json.dumps(dict(status='NOTHING_TO_DO', ready=ready)))
        return
    plan_path = ROOT/'00_manifest'/('sbds-plan-%s.json' % a.plan_tag)
    if plan_path.exists():
        raise SystemExit('refusing to overwrite plan: '+str(plan_path))
    subprocess.run([sys.executable, str(a.repo/'research/hprcv2/sbds_nrr/make_cohort_plan.py'),
        '--repo', str(a.repo), '--bindings', ROOT/'00_manifest/bindings',
        '--include', ','.join(pending),
        '--authorization', a.authorization or 'DRY-RUN-NOT-AUTHORIZED (driver default; not a launch)',
        '--scorer-dir', ROOT/'00_manifest/nrr-tooling-v1/scorer-v2',
        '--out', str(plan_path)], check=True)
    names_path = ROOT/'00_manifest'/('sbds-assemblies-%s.txt' % a.plan_tag)
    names_path.write_text(''.join(x+'\n' for x in pending))
    command = ['sbatch', '--array=0-%d%%%d' % (len(pending)-1, a.concurrency),
        '--cpus-per-task=1', '--mem=56G', '--time', a.time_limit,
        '--output', str(ROOT/'logs/sbds-%A_%a.log'),
        '--export=ALL,PLAN=%s,ASSEMBLIES_FILE=%s' % (plan_path, names_path),
        str(ROOT/'00_manifest/nrr-tooling-v1/submit_cohort.slurm')]
    print(json.dumps(dict(status='PREPARED', ready=len(ready), pending=len(pending),
        assemblies=pending, plan=str(plan_path), submit_command=' '.join(command))))
    if a.submit:
        if not a.authorization or 'DRY-RUN' in a.authorization:
            raise SystemExit('--submit requires a verbatim user authorization string')
        if a.concurrency > 32:
            raise SystemExit('concurrency above 32 requires new resource review')
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout)


if __name__ == '__main__':
    main()
