#!/usr/bin/env python3
"""Prepare an in-place verified binding for an assembly with trusted evidence.

For assemblies covered by the trusted index-recovery verification manifests,
creates an EXCLUSIVE stage directory containing symlinks to the original index
components plus a standard `crispr-binding.json` and a verified decompressed
reference. The normal per-run verification then hashes through the symlinks;
no original file is modified and no checksum is taken on trust.

The 454 assemblies WITHOUT trusted component evidence are rejected: preparing
them requires the full-content verification campaign first (documented in the
expanded-benchmark report). SBDS/PAM semantics are not involved here.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import time

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def verify_copy(source, dest, expected):
    if digest(source) != expected:
        raise ValueError('trusted checksum mismatch: '+str(source))
    with open(source, 'rb') as src, open(dest, 'xb') as dst:
        for block in iter(lambda: src.read(4*1024*1024), b''):
            dst.write(block)
    os.fsync(dst.fileno())
    if digest(dest) != expected:
        raise ValueError('copy verification failed: '+str(dest))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly', required=True)
    p.add_argument('--assemblies', required=True, help='trusted verified_assemblies.tsv')
    p.add_argument('--components', required=True, help='trusted verified_component_checksums.tsv')
    p.add_argument('--out', required=True, help='NEW exclusive stage directory')
    p.add_argument('--scoring-python', required=True, help='interpreter for reference_index.py')
    p.add_argument('--reference-index-helper', required=True, help='research/crispr_bulge/reference_index.py')
    a = p.parse_args()
    rows = [r for r in csv.DictReader(open(a.assemblies), delimiter='\t') if r['assembly_id'] == a.assembly]
    if len(rows) != 1:
        raise ValueError('assembly must appear exactly once in the trusted manifest')
    row = rows[0]
    components = {r['path']: r for r in csv.DictReader(open(a.components), delimiter='\t')
                  if r['assembly'] == a.assembly and r['status'] == 'READABLE'}
    out = Path(a.out)
    out.mkdir(mode=0o700)
    (out/'INCOMPLETE').write_text('Binding preparation not committed.\n')
    timings = {}
    start = time.monotonic()
    reference = out/'reference.fa'
    source_gz = Path(row['fasta'])
    h = hashlib.sha256()
    with gzip.open(source_gz, 'rb') as src, open(reference, 'xb') as dst:
        for block in iter(lambda: src.read(4*1024*1024), b''):
            dst.write(block)
            h.update(block)
        os.fsync(dst.fileno())
    reference_sha256 = h.hexdigest()
    if reference_sha256 != row['fasta_sha256']:
        raise ValueError('decompressed reference does not match trusted frozen FASTA hash')
    timings['decompress_seconds'] = time.monotonic()-start
    start = time.monotonic()
    for suffix in ('fai', 'gzi'):
        expected = row[f'{suffix}_sha256']
        if len(expected) != 64:
            raise ValueError('missing trusted '+suffix+' hash')
        verify_copy(str(source_gz)+'.'+suffix, str(reference)+'.'+suffix, expected)
    index_prefix = Path(row['index_prefix'])
    binding = dict(format='crispr-poc-index-binding-v1',
                   reference_sha256=reference_sha256, components_sha256={},
                   reference_check=dict(scope='trusted recovery verification + in-place symlink binding',
                                        bases=int(row['total_bases']), contigs=int(row['contigs'])))
    stage_prefix = out/'index'
    for suffix in SUFFIXES:
        original = str(index_prefix)+'.'+suffix
        entry = components.get(original)
        if entry is None or len(entry['sha256']) != 64:
            raise ValueError('no trusted component evidence for '+suffix)
        current = digest(original)
        if current != entry['sha256']:
            raise ValueError('original component changed since trusted scan: '+suffix)
        binding['components_sha256'][suffix] = current
        os.symlink(original, str(stage_prefix)+'.'+suffix)
    timings['component_verify_seconds'] = time.monotonic()-start
    with (str(stage_prefix)+'.crispr-binding.json').open('x') as f:
        json.dump(binding, f, indent=2, sort_keys=True)
        f.write('\n')
    start = time.monotonic()
    subprocess.run([a.scoring_python, '-I', a.reference_index_helper,
                    '--reference', str(reference), '--output', str(out/'reference-index.json')],
                   check=True, timeout=7200)
    timings['scoring_reference_index_seconds'] = time.monotonic()-start
    summary = dict(assembly=a.assembly, stage=str(out), reference=str(reference),
                   reference_sha256=reference_sha256, index=str(stage_prefix),
                   source_index=str(index_prefix), source_fasta=str(source_gz),
                   method='in-place symlink binding from trusted recovery evidence',
                   timings=timings,
                   limits='Full component checksums were re-verified against the trusted scan at bind time; per-run verification still hashes everything. Not an authenticity proof of the original build.')
    with (out/'binding.json').open('x') as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write('\n')
    (out/'INCOMPLETE').unlink()
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as e:
        print('binding preparation failed: '+str(e), file=sys.stderr)
        sys.exit(1)
