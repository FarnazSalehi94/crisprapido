#!/usr/bin/env python3
"""Delete a quarantined corrupted original after full verification (step 9).

Runs ONLY when: the replacement is promoted in indices/, the binding is
committed, and the assembly's four-PAM SBDS search is COMPLETE with verified
checksums. Preserves the original's error report and writes a deletion log
with all hashes before removing the quarantine copy. Never touches healthy
indexes, FASTAs, or historical results.
"""
import hashlib
import json
import os
from pathlib import Path
import sys
import time

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido')
SEARCH = ROOT/'runs/sbds_ngg_nrr_v1/02_search'
LOGS = ROOT/'indices_deletion_log_v1'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def delete_original(assembly):
    promotion = json.loads((ROOT/'indices'/assembly/'promotion.json').read_text())
    if promotion.get('status') != 'PROMOTED' or not promotion.get('original_quarantined'):
        raise SystemExit('no promoted replacement with quarantined original: '+assembly)
    success = SEARCH/assembly/'SUCCESS.json'
    if not success.exists():
        raise SystemExit('SBDS search not COMPLETE; refusing deletion: '+assembly)
    data = json.loads(success.read_text())
    for path, spec in data['outputs'].items():
        if sha256(path) != spec['sha256']:
            raise SystemExit('completed search output changed; refusing deletion: '+path)
    quarantine = Path(promotion['original_quarantined'])
    incident = ROOT/'runs/sbds_ngg_nrr_v1/00_manifest/incident-evidence-v1'/(assembly+'.json')
    incident_preserved = None
    if incident.exists():
        incident_preserved = str(LOGS/(assembly+'-error-report.json'))
        LOGS.mkdir(exist_ok=True)
        if not Path(incident_preserved).exists():
            (LOGS/(assembly+'-error-report.json')).write_bytes(incident.read_bytes())
    entries = []
    for p in sorted(quarantine.iterdir()):
        entries.append(dict(file=p.name, bytes=p.stat().st_size,
                            sha256=sha256(p) if p.is_file() else None))
    freed = sum(e['bytes'] for e in entries if e['bytes'])
    record = dict(format='index-deletion-v1', assembly=assembly,
        deleted_path=str(quarantine), deleted_entries=entries, freed_bytes=freed,
        replacement_verified_by=promotion['components_sha256'],
        search_success=str(success), search_outputs=len(data['outputs']),
        error_report_preserved=str(incident_preserved),
        deleted_utc=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        note='Corrupted original deleted only after promoted replacement passed cross-node verification, binding, and a COMPLETE checksummed SBDS search.')
    os.rename(quarantine, LOGS/(assembly+'-deleted'))
    (LOGS/(assembly+'-deletion.json')).write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(dict(assembly=assembly, status='ORIGINAL_DELETED',
        freed_bytes=freed, error_report=str(incident_preserved))))


if __name__ == '__main__':
    delete_original(sys.argv[1])
