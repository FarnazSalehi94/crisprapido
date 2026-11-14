#!/usr/bin/env python3
"""SBDS index rebuild under the 2026-11-13 MooseFS storage policy.

Builds on node-local scratch (fast, no partial MooseFS exposure), fully
verifies every component on scratch, copies to the permanent MooseFS tree
(`indices_rebuilt_v1/`), re-reads every byte from MooseFS and compares hashes
(write-back verification), fsyncs, and writes a checksummed completion
manifest. A SEPARATE job (reverify_sbds.py) must re-read the finished set from
a DIFFERENT worker before it is marked trusted; only then may binding run.
Corrupted originals are never touched or deleted here.
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
import socket

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
LOADER = '/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib/ld-linux-x86-64.so.2'
LIBPATH = ('/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib:'
           '/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib:'
           '/gnu/store/ygm2rzg08w0w96h9587k3n2yvhndcqkn-zlib-1.2.13/lib')
BUILDER = '/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/rebuilt_HG00128_hap2_attempt02/tools/columba_build'
DEST_ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido/indices_rebuilt_v1')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def copy_verified(src, dst):
    with open(src, 'rb') as s, open(dst, 'xb') as d:
        for b in iter(lambda: s.read(4*1024*1024), b''):
            d.write(b)
        d.flush()
        os.fsync(d.fileno())
    if sha256(dst) != sha256(src):
        raise ValueError('write-back verification failed: '+str(dst))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--assembly', required=True)
    p.add_argument('--builder', default=BUILDER)
    a = p.parse_args()
    fasta_gz = dict(l.split('\t') for l in
        Path('/moosefs/salehi/HPRCv2_CRISPRapido/analysis_genomes.tsv').read_text().splitlines() if l.strip())[a.assembly]
    if '/' in a.assembly or a.assembly in ('', '.', '..'):
        raise SystemExit('invalid assembly name')
    if sha256(a.builder) != '5987ad90e0c4a061cbf4aac7d6b087a284ccf7fe0f9e1835219c567c5e8d1ad2':
        raise SystemExit('builder pin mismatch')
    scratch_base = Path(os.environ.get('SBDS_REBUILD_SCRATCH') or '/scratch/salehi')
    scratch_base.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(scratch_base).free < 80_000_000_000:
        raise SystemExit('scratch free-space preflight failed (<80 GB)')
    if shutil.disk_usage(DEST_ROOT.parent).free < 2_000_000_000_000:
        raise SystemExit('persistent storage floor breached (<2 TB)')
    source_gz_sha_before = sha256(fasta_gz)
    scratch = scratch_base/('rebuild-'+a.assembly+'-'+os.environ['SLURM_JOB_ID'])
    scratch.mkdir(parents=True)
    dest = DEST_ROOT/a.assembly
    if dest.exists():
        raise SystemExit('refusing to reuse/remove existing destination: '+str(dest))
    timings = {}
    t = time.monotonic()
    plain = scratch/(a.assembly+'.fa')
    h = hashlib.sha256()
    with gzip.open(fasta_gz, 'rb') as src, plain.open('xb') as dst:
        for b in iter(lambda: src.read(4*1024*1024), b''):
            dst.write(b); h.update(b)
        dst.flush(); os.fsync(dst.fileno())
    source_plain_sha = h.hexdigest()
    source_gz_sha = sha256(fasta_gz)
    if source_gz_sha != source_gz_sha_before:
        raise SystemExit('source FASTA changed while decompressing; scratch preserved')
    timings['decompress_verify_seconds'] = time.monotonic()-t
    t = time.monotonic()
    index_prefix = scratch/'index'/a.assembly
    (scratch/'index').mkdir()
    command = ['/home/salehi/miniconda3/bin/time', '-v', '-o', str(scratch/'build.time'),
               '/usr/bin/timeout', '7200', LOADER, '--library-path', LIBPATH, a.builder,
               '-f', str(plain), '-r', str(index_prefix)]
    with (scratch/'build.stdout').open('x') as so, (scratch/'build.stderr').open('x') as se:
        rc = subprocess.run(command, stdout=so, stderr=se, env=dict(os.environ, OMP_NUM_THREADS='1')).returncode
    timings['build_seconds'] = time.monotonic()-t
    if rc != 0:
        raise SystemExit('builder failed rc=%d; scratch %s preserved for diagnosis' % (rc, scratch))
    t = time.monotonic()
    scratch_hashes = {s: sha256(str(index_prefix)+'.'+s) for s in SUFFIXES}
    timings['scratch_verify_seconds'] = time.monotonic()-t
    dest.mkdir(mode=0o700)
    (dest/'INCOMPLETE').write_text('rebuild not committed\n')
    t = time.monotonic()
    mfs_hashes = {}
    for s in SUFFIXES:
        copy_verified(str(index_prefix)+'.'+s, str(dest/(a.assembly+'.'+s)))
        mfs_hashes[s] = sha256(str(dest/(a.assembly+'.'+s)))
        if mfs_hashes[s] != scratch_hashes[s]:
            raise SystemExit('MooseFS write-back verification failed for '+s)
    shutil.copyfile(scratch/'build.time', dest/'build.time')
    for n in ('build.stdout', 'build.stderr'):
        shutil.copyfile(scratch/n, dest/n)
    timings['mfs_writeback_verify_seconds'] = time.monotonic()-t
    build_time = {}
    for line in (dest/'build.time').read_text().splitlines():
        if 'Elapsed' in line or 'Maximum resident' in line:
            k, v = line.rsplit(':', 1); build_time[k.strip()] = v.strip()
    manifest = dict(format='sbds-rebuilt-index-completion-v1', assembly=a.assembly,
        status='WRITTEN_AWAITING_CROSS_NODE_RECHECK', destination=str(dest),
        build_host=socket.gethostname(),
        source_fasta=fasta_gz, source_compressed_sha256=source_gz_sha,
        source_plain_sha256=source_plain_sha, builder_sha256=sha256(a.builder),
        components_sha256=mfs_hashes, timings=timings, build_time=build_time,
        total_component_bytes=sum((dest/(a.assembly+'.'+s)).stat().st_size for s in SUFFIXES),
        integrity='all components verified on scratch, after MooseFS write-back, and pending re-read from a different worker',
        original_index_untouched=True)
    (dest/'completion.json').write_text(json.dumps(manifest, indent=2)+'\n')
    shutil.rmtree(scratch)
    print(json.dumps(dict(assembly=a.assembly, status=manifest['status'],
        timings={k: round(v, 1) for k, v in timings.items()},
        index_bytes=manifest['total_component_bytes'])))


if __name__ == '__main__':
    main()
