#!/usr/bin/env python3
"""Promote a TRUSTED rebuilt index into the authoritative indices/ directory.

Safe, recoverable directory switch (same filesystem renames only):
  1. staging must be TRUSTED_CROSS_NODE with a completion manifest
  2. rename indices/<asm> -> indices_quarantine_v1/<asm>-original (rollback point)
  3. rename indices_rebuilt_v1/<asm> -> indices/<asm>
  4. write promotion.json (rollback instructions, hashes)
Any failure rolls both renames back and stops that replacement. The
quarantined original is preserved until the assembly's SBDS search completes
and delete_original.py runs with its error report and deletion log.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido')
INDICES = ROOT/'indices'
STAGING = ROOT/'indices_rebuilt_v1'
QUARANTINE = ROOT/'indices_quarantine_v1'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def promote(assembly):
    staging = STAGING/assembly
    dest = INDICES/assembly
    manifest = json.loads((staging/'completion.json').read_text())
    if manifest['status'] != 'TRUSTED_CROSS_NODE':
        raise SystemExit('staging set not cross-node trusted: '+assembly)
    quarantine = QUARANTINE/(assembly+'-original')
    if quarantine.exists():
        raise SystemExit('quarantine entry already exists: '+str(quarantine))
    QUARANTINE.mkdir(exist_ok=True)
    had_original = dest.exists()
    promoted = False
    try:
        if had_original:
            os.rename(dest, quarantine)
        os.rename(staging, dest)
        promoted = True
        # verify the promoted set in place against the completion manifest
        for s in SUFFIXES:
            if sha256(dest/(assembly+'.'+s)) != manifest['components_sha256'][s]:
                raise SystemExit('post-promotion verification failed: '+s)
        record = dict(format='index-promotion-v1', assembly=assembly,
            status='PROMOTED', promoted_from=str(staging), authoritative=str(dest),
            original_quarantined=str(quarantine) if had_original else None,
            original_existed=had_original,
            rollback=('mv %s %s' % (dest, staging)) + (' && mv %s %s' % (quarantine, dest) if had_original else ''),
            components_sha256=manifest['components_sha256'],
            deletion_precondition='SBDS search SUCCESS for this assembly plus reviewed deletion log',
            note='original preserved in quarantine; delete only via delete_original.py')
        (dest/'promotion.json').write_text(json.dumps(record, indent=2)+'\n')
        print(json.dumps(dict(assembly=assembly, status='PROMOTED',
            original_quarantined=had_original)))
    except BaseException:
        if promoted:
            os.rename(dest, staging)
        if had_original and quarantine.exists() and not dest.exists():
            os.rename(quarantine, dest)
        raise


if __name__ == '__main__':
    promote(sys.argv[1])
