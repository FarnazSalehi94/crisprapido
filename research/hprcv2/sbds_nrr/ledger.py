#!/usr/bin/env python3
"""Audit trusted reference/index-binding evidence for the full 466-input cohort.

Metadata-only cross-check of existing manifests; NO index bytes are read and
no verdict here authorizes a search. Deep full-content verification is a
separate, explicitly budgeted campaign (prepare_binding --deep).
"""
import argparse
import csv
import json
from pathlib import Path

RECOVERY = Path('/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/recovery_report')
SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cohort', type=Path, default=Path('/moosefs/salehi/HPRCv2_CRISPRapido/analysis_genomes.tsv'))
    p.add_argument('--recovery', type=Path, default=RECOVERY)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    cohort = [line.split('\t')[0] for line in a.cohort.read_text().splitlines() if line.strip()]
    if len(cohort) != len(set(cohort)) or len(cohort) != 466:
        raise ValueError('cohort manifest must list 466 unique assemblies')
    verified = {r['assembly_id']: r for r in csv.DictReader(open(a.recovery/'verified_assemblies.tsv'), delimiter='\t')}
    components = {}
    for r in csv.DictReader(open(a.recovery/'verified_component_checksums.tsv'), delimiter='\t'):
        if r['status'] == 'READABLE' and r['matches_readable_original'] != 'False':
            components.setdefault(r['assembly'], {})[r['path']] = r['sha256']
    rows = []
    for assembly in cohort:
        v = verified.get(assembly)
        if v is None:
            rows.append(dict(assembly=assembly, ledger_status='UNVERIFIED',
                reason='no trusted full-content index evidence; deep verification required'))
            continue
        missing = [s for s in SUFFIXES
                   if (str(v['index_prefix'])+'.'+s) not in components.get(assembly, {})
                   and v['recovery_method'] != 'ORIGINAL_VERIFIED']
        ready = Path(str(v['index_prefix'])+'.bwt').exists()
        rows.append(dict(assembly=assembly, ledger_status=v['verification_status'],
            method=v['recovery_method'], index_prefix=v['index_prefix'],
            frozen_fasta=v['frozen_fasta'], component_evidence='full checksums' if not missing else 'incomplete:'+','.join(missing),
            index_present=ready,
            reason='' if ready else 'verified index prefix missing on disk'))
    bound = {r['assembly'] for r in rows if r['ledger_status'] in ('ORIGINAL_VERIFIED', 'REBUILT_VERIFIED')}
    committed = {}
    runs = Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1')
    bindings = runs/'00_manifest/bindings'
    batch = runs/'00_manifest/verification-batch-1'
    for stage in sorted(bindings.glob('*/binding.json')) if bindings.exists() else []:
        committed[stage.parent.name] = 'DEEP_VERIFIED_BOUND'
    batch_failures = {}
    if batch.exists():
        for rec in sorted(batch.glob('*.json')):
            data = json.loads(rec.read_text())
            if data['status'] == 'FAILED':
                batch_failures[data['assembly']] = data['failure_kind']+': '+data['error'][:80]
    report = dict(status='METADATA_AUDIT', cohort=466,
        verified_bindings=sorted(bound), verified_count=len(bound),
        unverified_count=466-len(bound),
        unverified=sorted(r['assembly'] for r in rows if r['ledger_status'] == 'UNVERIFIED'),
        committed_bindings=committed, batch_failures=batch_failures,
        ready_binding_sidecars={
            'HG00735_pat': '/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1',
            'HG01496_pat': '/export/local/home/salehi/crisprapido-hprc-import-hg01496-20261112-v1'},
        rows=rows,
        caveats=['Metadata presence is never content verification; deep mode hashes every component byte.',
                 'Verified sets were confirmed 2026-10-28/29; components could degrade again before use.',
                 'This audit authorizes nothing; per-run binding checks still re-verify at execution time.'])
    a.out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(verified=len(bound), unverified=466-len(bound))))

if __name__ == '__main__':
    main()
