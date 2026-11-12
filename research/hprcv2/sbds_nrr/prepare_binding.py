#!/usr/bin/env python3
"""SBDS binding preparation and full-content verification; never a search.

Repairs the withdrawn Brunello draft: correct open files, correct digest
semantics (ledger fasta_sha256 is the COMPRESSED .gz hash), fsync before close,
no string .open(), exclusive stage, and a --deep mode for assemblies without
trusted evidence (full component reads produce new checksums; that is the
only place multi-TB reads happen, one assembly at a time, bounded by the
caller's Slurm budget). Reference decompression and scoring index are
per-assembly one-time costs.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
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
        dst.flush()
        os.fsync(dst.fileno())
    if digest(dest) != expected:
        raise ValueError('copy verification failed: '+str(dest))


def decompress_reference(source_gz, reference, expected_gz_sha256):
    h_gz = hashlib.sha256()
    with open(source_gz, 'rb') as raw:
        for block in iter(lambda: raw.read(4*1024*1024), b''):
            h_gz.update(block)
    if expected_gz_sha256 and h_gz.hexdigest() != expected_gz_sha256:
        raise ValueError('frozen compressed FASTA changed since trusted scan')
    h = hashlib.sha256()
    with gzip.open(source_gz, 'rb') as src, open(reference, 'xb') as dst:
        for block in iter(lambda: src.read(4*1024*1024), b''):
            dst.write(block)
            h.update(block)
        dst.flush()
        os.fsync(dst.fileno())
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly', required=True)
    p.add_argument('--index-prefix', required=True)
    p.add_argument('--fasta-gz', required=True)
    p.add_argument('--fai', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--scoring-python', required=True)
    p.add_argument('--reference-index-helper', required=True)
    p.add_argument('--trusted', action='store_true', help='require trusted component checksums')
    p.add_argument('--component-manifest', default='', help='trusted verified_component_checksums.tsv for drift comparison')
    p.add_argument('--fasta-gz-sha256', default='', help='trusted compressed FASTA hash')
    p.add_argument('--fai-sha256', default='', help='trusted fai hash')
    a = p.parse_args()
    out = Path(a.out)
    if out.exists():
        if (out/'binding.json').exists():
            print(json.dumps({'status': 'ALREADY_BOUND', 'assembly': a.assembly, 'stage': str(out)}))
            return
        raise ValueError('refusing to reuse/remove an existing stage directory: '+str(out))
    out.mkdir(mode=0o700)
    (out/'INCOMPLETE').write_text('Binding preparation not committed.\n')
    timings = {}
    start = time.monotonic()
    reference = out/'reference.fa'
    reference_sha256 = decompress_reference(a.fasta_gz, reference, a.fasta_gz_sha256)
    timings['decompress_verify_seconds'] = time.monotonic()-start
    start = time.monotonic()
    if a.fai_sha256:
        verify_copy(a.fai, str(reference)+'.fai', a.fai_sha256)
    else:
        shutil.copyfile(a.fai, str(reference)+'.fai')
        (out/'fai-unverified.txt').write_text('No trusted fai hash supplied; fai copied but content-unverified.\n')
    stage_prefix = out/'index'
    components = {}
    trusted = {}
    if a.component_manifest:
        for r in csv.DictReader(open(a.component_manifest), delimiter='\t'):
            if r['status'] == 'READABLE' and r['matches_readable_original'] != 'False':
                trusted[r['path']] = r['sha256']
    start = time.monotonic()
    for suffix in SUFFIXES:
        original = str(Path(a.index_prefix))+'.'+suffix
        h = hashlib.sha256()
        with open(original, 'rb') as f:
            for block in iter(lambda: f.read(4*1024*1024), b''):
                h.update(block)
        components[suffix] = h.hexdigest()
        if a.trusted:
            expected = trusted.get(original)
            if expected is None:
                raise ValueError('no trusted component evidence for '+suffix)
            if expected != components[suffix]:
                raise ValueError('component drifted since trusted scan: '+suffix)
        os.symlink(original, str(stage_prefix)+'.'+suffix)
    timings['component_full_read_seconds'] = time.monotonic()-start
    start = time.monotonic()
    with open(str(stage_prefix)+'.crispr-binding.json', 'x') as f:
        json.dump(dict(format='crispr-poc-index-binding-v1', reference_sha256=reference_sha256,
                       components_sha256=components,
                       reference_check=dict(scope=('trusted-manifest re-verification' if a.trusted else
                                                   'first full-content verification; new checksums produced here'),
                                            verification='all component bytes read and hashed this run')),
                  f, indent=2, sort_keys=True)
        f.write('\n')
    start = time.monotonic()
    env = dict(os.environ, PYTHONPATH=str(Path(a.reference_index_helper).parent))
    subprocess.run([a.scoring_python, a.reference_index_helper,
                    '--reference', str(reference), '--output', str(out/'reference-index.json')],
                   check=True, timeout=7200, env=env)
    timings['scoring_reference_index_seconds'] = time.monotonic()-start
    with (out/'binding.json').open('x') as f:
        json.dump(dict(assembly=a.assembly, stage=str(out), reference=str(reference),
                       reference_sha256=reference_sha256, index=str(stage_prefix),
                       source_index=a.index_prefix, source_fasta=a.fasta_gz,
                       trusted=a.trusted, timings=timings,
                       limits='Binding re-verification still runs per search job; this is not an authenticity proof of the original build.'),
                  f, indent=2, sort_keys=True)
        f.write('\n')
    (out/'INCOMPLETE').unlink()
    print(json.dumps(dict(status='BOUND', assembly=a.assembly, stage=str(out), timings=timings)))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as e:
        print('binding preparation failed: '+str(e), file=sys.stderr)
        sys.exit(1)
