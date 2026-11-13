#!/usr/bin/env python3
"""Migrate the three /export pilot rebuilt indexes into the MooseFS tree.

Writes completion manifests and performs full write-back verification; the
separate reverify_sbds.py step (different worker) marks them trusted. The
/export originals are preserved until the migrated copies pass cross-node
recheck and a deletion review; nothing is deleted by this script.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
SRC = Path('/export/local/home/salehi/sbds_rebuilt_indexes_v1')
DEST_ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido/indices_rebuilt_v1')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def migrate(assembly):
    src = SRC/assembly/'index'/assembly
    rebuild_meta = json.loads((SRC/assembly/'rebuild.json').read_text())
    dest = DEST_ROOT/assembly
    if dest.exists():
        raise SystemExit('destination exists: '+str(dest))
    dest.mkdir(mode=0o700)
    (dest/'INCOMPLETE').write_text('migration not committed\n')
    t = time.monotonic()
    hashes = {}
    for s in SUFFIXES:
        out = dest/(assembly+'.'+s)
        with open(str(src)+'.'+s, 'rb') as i, open(out, 'xb') as o:
            for b in iter(lambda: i.read(4*1024*1024), b''):
                o.write(b)
            o.flush(); os.fsync(o.fileno())
        if sha256(out) != sha256(str(src)+'.'+s):
            raise SystemExit('write-back verification failed: '+s)
        hashes[s] = sha256(out)
    if hashes != rebuild_meta['components_sha256']:
        raise SystemExit('migrated hashes disagree with pilot manifest')
    for n in ('build.time', 'build.stdout', 'build.stderr'):
        shutil.copyfile(SRC/assembly/n, dest/n)
    manifest = dict(format='sbds-rebuilt-index-completion-v1', assembly=assembly,
        status='WRITTEN_AWAITING_CROSS_NODE_RECHECK', destination=str(dest),
        source_fasta=rebuild_meta['source_fasta'],
        source_compressed_sha256=rebuild_meta['source_compressed_sha256'],
        source_plain_sha256=rebuild_meta['source_plain_sha256'],
        builder_sha256=rebuild_meta['builder_sha256'],
        components_sha256=hashes, timings=rebuild_meta['timings'],
        build_time={},
        total_component_bytes=sum((dest/(assembly+'.'+s)).stat().st_size for s in SUFFIXES),
        provenance='migrated from /export pilot (jobs 2871822); hashes verified identical',
        original_index_untouched=True)
    (dest/'completion.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return dict(assembly=assembly, seconds=round(time.monotonic()-t, 1),
        index_bytes=manifest['total_component_bytes'])


if __name__ == '__main__':
    results = [migrate(a) for a in ('HG005_pat', 'HG03225_hap1', 'HG03784_hap1')]
    print(json.dumps(results))
