#!/usr/bin/env python3
"""Re-point an existing committed binding at the promoted authoritative index.

Preserves the binding sidecar byte-for-byte (its content is location-
independent: component hashes + reference hash), so completed SBDS searches
keep their valid scientific signatures and are never rerun. Verifies every
component at the new location against the binding record before switching.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def repoint(assembly):
    binding_dir = ROOT/'runs/sbds_ngg_nrr_v1/00_manifest/bindings'/assembly
    binding = json.loads((binding_dir/'index.crispr-binding.json').read_text())
    before = sha256(binding_dir/'index.crispr-binding.json')
    dest_prefix = ROOT/'indices'/assembly/assembly
    for s in SUFFIXES:
        new_path = str(dest_prefix)+'.'+s
        if sha256(new_path) != binding['components_sha256'][s]:
            raise SystemExit('component mismatch at authoritative location: '+s)
    for s in SUFFIXES:
        link = binding_dir/('index.'+s)
        new_path = str(dest_prefix)+'.'+s
        if link.is_symlink() and os.readlink(link) == new_path:
            continue
        tmp = binding_dir/('index.'+s+'.new')
        if tmp.is_symlink() or tmp.exists():
            tmp.unlink()
        os.symlink(new_path, tmp)
        os.replace(tmp, link)
    if sha256(binding_dir/'index.crispr-binding.json') != before:
        raise SystemExit('binding sidecar changed unexpectedly')
    (binding_dir/'repoint-record.json').write_text(json.dumps(dict(
        assembly=assembly, authoritative_index_prefix=str(dest_prefix),
        sidecar_unchanged=True, note='symlinks re-pointed to promoted authoritative index'), indent=2)+'\n')
    print(json.dumps(dict(assembly=assembly, status='REPOINTED', sidecar_unchanged=True)))


if __name__ == '__main__':
    repoint(sys.argv[1])
