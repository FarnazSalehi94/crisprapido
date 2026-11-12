#!/usr/bin/env python3
"""Controlled single-assembly index REBUILD pilot onto healthy local storage.

Rebuilds one failed assembly's Columba index from its verified frozen FASTA
onto /export (local RAID), never onto the degraded MooseFS index store, and
never touching the original index. Requires the incident scan's verified
source-FASTA evidence for the assembly. Records builder wall time, GNU-time
peak RSS, storage, and full component verification. The rebuilt index is then
eligible for binding via the standard deep verification helper.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1')
LOADER = '/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib/ld-linux-x86-64.so.2'
LIBPATH = ('/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib:'
           '/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib:'
           '/gnu/store/ygm2rzg08w0w96h9587k3n2yvhndcqkn-zlib-1.2.13/lib')
BUILDER = '/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/rebuilt_HG00128_hap2_attempt02/tools/columba_build'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--assembly', required=True)
    p.add_argument('--out', type=Path, required=True, help='NEW exclusive /export stage')
    a = p.parse_args()
    incident = json.loads((ROOT/'00_manifest/incident-evidence-v1'/(a.assembly+'.json')).read_text())
    if incident['source_fasta']['status'] != 'READABLE':
        raise SystemExit('refusing rebuild: source FASTA unreadable for '+a.assembly)
    fasta_gz = incident['source_fasta']['path']
    out = a.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise SystemExit('refusing to reuse/remove existing stage: '+str(out))
    out.mkdir(mode=0o700)
    (out/'INCOMPLETE').write_text('rebuild pilot not committed\n')
    started = time.monotonic()
    (out/'input').mkdir(exist_ok=True)
    plain = out/'input'/(a.assembly+'.fa')
    h = hashlib.sha256()
    with gzip.open(fasta_gz, 'rb') as src, plain.open('xb') as dst:
        for b in iter(lambda: src.read(4*1024*1024), b''):
            dst.write(b); h.update(b)
        dst.flush(); os.fsync(dst.fileno())
    if h.hexdigest() != incident['source_fasta']['decompressed_sha256']:
        raise SystemExit('decompressed source hash drift vs incident evidence')
    decompress_seconds = time.monotonic()-started
    command = ['/home/salehi/miniconda3/bin/time', '-v', '-o', str(out/'build.time'),
               '/usr/bin/timeout', '7200', LOADER, '--library-path', LIBPATH, BUILDER,
               '-f', str(plain), '-r', str(out/'index'/a.assembly)]
    (out/'index').mkdir()
    build_started = time.monotonic()
    with (out/'build.stdout').open('x') as so, (out/'build.stderr').open('x') as se:
        rc = subprocess.run(command, stdout=so, stderr=se, env=dict(os.environ, OMP_NUM_THREADS='1')).returncode
    build_seconds = time.monotonic()-build_started
    metrics = dict(wall_seconds=build_seconds, returncode=rc,
                   command=command[:2]+['...']+[BUILDER, '-f', plain.name, '-r', 'index/'+a.assembly])
    if rc != 0:
        raise SystemExit('builder failed rc=%d; stage preserved at %s' % (rc, out))
    # Full verification of the rebuilt set: every component readable + hashed.
    verify_started = time.monotonic()
    components = {}
    for s in SUFFIXES:
        components[s] = sha256(str(out/'index'/a.assembly)+'.'+s)
    verify_seconds = time.monotonic()-verify_started
    total_bytes = sum((out/'index'/(a.assembly+'.'+s)).stat().st_size for s in SUFFIXES)
    report = dict(assembly=a.assembly, status='REBUILT_VERIFIED',
        stage=str(out), source_fasta=fasta_gz,
        source_compressed_sha256=sha256(fasta_gz),
        source_plain_sha256=incident['source_fasta']['decompressed_sha256'],
        builder_sha256=sha256(BUILDER), loader=LOADER,
        original_index_untouched='/moosefs/salehi/HPRCv2_CRISPRapido/indices/'+a.assembly,
        timings=dict(decompress_seconds=decompress_seconds, build_seconds=build_seconds,
                     component_verify_seconds=verify_seconds),
        index_component_bytes=total_bytes, components_sha256=components,
        note='Rebuilt onto local /export RAID from verified frozen FASTA; NOT onto the degraded MooseFS index store.')
    (out/'rebuild.json').write_text(json.dumps(report, indent=2)+'\n')
    (out/'INCOMPLETE').unlink()
    print(json.dumps({k: report[k] for k in ('assembly', 'status', 'timings', 'index_component_bytes')}))


if __name__ == '__main__':
    main()
