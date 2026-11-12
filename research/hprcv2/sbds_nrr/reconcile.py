#!/usr/bin/env python3
"""Reconcile all verification evidence into ONE authoritative cohort manifest.

Exactly 466 unique assemblies with mutually exclusive statuses:
  READY_BOUND            - committed binding.json (deep- or trusted-verified)
  TRUSTED_EVIDENCE       - trusted recovery evidence, no committed binding yet
  INDEX_EIO              - deep verification failed on unreadable components
  UNATTEMPTED            - no deep verification attempt yet
The earlier '20 search-ready' figure double-counted HG02055_mat and
NA18906_mat (both trusted AND bound); unique ready = 18.
"""
import argparse
import json
from pathlib import Path

RECOVERY = Path('/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/recovery_report')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cohort', type=Path, default=Path('/moosefs/salehi/HPRCv2_CRISPRapido/analysis_genomes.tsv'))
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    cohort = [l.split('\t')[0] for l in a.cohort.read_text().splitlines() if l.strip()]
    if len(cohort) != 466 or len(set(cohort)) != 466:
        raise ValueError('cohort must contain 466 unique assemblies, found %d rows / %d unique'
                         % (len(cohort), len(set(cohort))))
    runs = Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/00_manifest')
    import csv
    trusted = {r['assembly_id'] for r in csv.DictReader(open(RECOVERY/'verified_assemblies.tsv'), delimiter='\t')}
    bound = {s.parent.name for s in (runs/'bindings').glob('*/binding.json') if not (s.parent/'INCOMPLETE').exists()}
    eio = {}
    batch = runs/'verification-batch-1'
    if batch.exists():
        for rec in sorted(batch.glob('*.json')):
            data = json.loads(rec.read_text())
            if data['status'] == 'FAILED':
                eio[data['assembly']] = data['failure_kind']
    rows = []
    for assembly in cohort:
        if assembly in bound:
            status = 'READY_BOUND'
        elif assembly in trusted:
            status = 'TRUSTED_EVIDENCE'
        elif assembly in eio:
            status = 'INDEX_EIO' if eio[assembly] == 'EIO' else 'INDEX_FAILED_'+eio[assembly]
        else:
            status = 'UNATTEMPTED'
        rows.append(dict(assembly=assembly, status=status))
    from collections import Counter
    counts = Counter(r['status'] for r in rows)
    assert sum(counts.values()) == 466 and len(rows) == len(set(r['assembly'] for r in rows))
    report = dict(format='sbds-cohort-manifest-v1', cohort_size=466,
        status_counts=dict(counts), mutually_exclusive=True,
        reconciliation=dict(
            previous_reported_ready=20,
            double_counted=sorted(bound & trusted),
            corrected_ready=len(bound | trusted),
            explanation='READY_BOUND and TRUSTED_EVIDENCE are disjoint; every assembly appears exactly once.'),
        rows=rows,
        caveats=['READY_BOUND requires a committed binding; TRUSTED_EVIDENCE assemblies still need binding preparation.',
                 'INDEX_EIO originals are preserved untouched; status refers to index components only, never the assembly FASTA.',
                 'UNATTEMPTED means no full-content evidence either way.'])
    a.out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(counts=dict(counts), corrected_ready=len(bound | trusted), double_counted=sorted(bound & trusted))))


if __name__ == '__main__':
    main()
