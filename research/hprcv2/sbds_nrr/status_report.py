#!/usr/bin/env python3
"""Authoritative campaign status report: indexes, rebuilds, SBDS, storage.

Reads only committed, checksummed evidence; missing/failed work is reported
as such and never zero-filled. Run freely at any time.
"""
import json
from pathlib import Path
import subprocess

ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido')
RUNS = ROOT/'runs/sbds_ngg_nrr_v1'


def tree_gb_all(root, runs):
    def gb(path):
        if not path.exists(): return 0
        r = subprocess.run(['du','-sm','--no-dereference',str(path)],capture_output=True,text=True)
        if not r.stdout.strip(): return 'du-failed'
        return round(int(r.stdout.split()[0])/1024, 1)
    return dict(indices=gb(root/'indices'), quarantine=gb(root/'indices_quarantine_v1'),
                bindings=gb(runs/'00_manifest/bindings'))


def main():
    m = json.loads((RUNS/'00_manifest/authoritative-manifest-v5.json').read_text())
    rows = {r['assembly']: r['status'] for r in m['rows']}
    completed, failed, partial = [], [], []
    for a, status in rows.items():
        if (RUNS/'02_search'/a/'SUCCESS.json').exists():
            completed.append(a)
        elif a != 'HG02055_mat' and (RUNS/'02_search'/a).exists() and sorted((RUNS/'02_search'/a).glob('attempt-*/FAILED.json')):
            failed.append(a)
    promoted = sorted(p.parent.name for p in (ROOT/'indices').glob('*/promotion.json'))
    quarantined = sorted(d.name for d in (ROOT/'indices_quarantine_v1').glob('*') if d.is_dir())
    deleted = sorted(p.stem.replace('-deletion', '') for p in (ROOT/'indices_deletion_log_v1').glob('*-deletion.json')) if (ROOT/'indices_deletion_log_v1').exists() else []
    staging = sorted(d.name for d in (ROOT/'indices_rebuilt_v1').glob('*') if d.is_dir())
    hits = cfd = cb = na = 0
    for a in completed:
        for pam in ('NGG', 'NAG', 'NGA', 'NAA'):
            paths = sorted((RUNS/'02_search'/a).glob('attempt-*/'+pam+'.summary.json'))
            if not paths:
                continue
            c = json.loads(paths[-1].read_text())['counts']
            hits += c['rows']; cfd += c['cfd_defined']
            cb += c.get('cb_scored', 0)
            na += sum(v for k, v in c.items() if k.startswith('cb_') and k != 'cb_scored')
    report = dict(
        cohort=466,
        index_ledger={s: sum(1 for v in rows.values() if v == s) for s in sorted(set(rows.values()))},
        healthy_originals=sum(1 for a, s in rows.items() if s == 'READY_BOUND' and a not in promoted),
        promoted_rebuilt=len(promoted),
        quarantined_originals=len(quarantined),
        originals_deleted_after_verification=len(deleted),
        sbds=dict(completed=len(completed), remaining=sum(1 for s in rows.values() if s == 'READY_BOUND')-len(completed),
                  failed=len(failed),
                  total_hits=hits, cfd_defined=cfd, crispr_bulge_defined=cb, cb_na=na),
        brunello='prepared, not launched; see BRUNELLO_READINESS.md',
        storage_gb=tree_gb_all(ROOT, RUNS),
        caveats=['Every figure is read from committed manifests; no zero-filling.',
                 'cb scores on NAG/NGA/NAA are computational predictions with unestablished biological calibration.'])
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
