#!/usr/bin/env python3
"""Recover only the manifest of build 2875899; never rebuild or delete evidence."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import socket
from rebuild_sbds import SUFFIXES, BUILDER, sha256


def verify(source, scratch, dest, expected):
    if sha256(source) != expected:
        raise ValueError('derived source digest mismatch')
    h = hashlib.sha256()
    with gzip.open(source, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''):
            h.update(b)
    if sha256(scratch/'GRCh38.fa') != h.hexdigest():
        raise ValueError('scratch reference differs from derived source')
    hashes = {}
    for suffix in SUFFIXES:
        name = 'GRCh38.'+suffix
        hashes[suffix] = sha256(dest/name)
        if hashes[suffix] != sha256(scratch/'index'/name):
            raise ValueError('scratch/staged component mismatch: '+suffix)
    return h.hexdigest(), hashes


def main():
    root = Path('/moosefs/salehi/HPRCv2_CRISPRapido')
    derived = root/'runs/sbds_ngg_nrr_v1/00_manifest/grch38-acgtn-fix-v2'
    dest = root/'indices_rebuilt_v1/GRCh38'
    scratch = Path('/scratch/salehi/rebuild-GRCh38-2875899')
    if socket.gethostname().split('.')[0] != 'octopus10':
        raise ValueError('must verify original build scratch on octopus10')
    if (dest/'completion.json').exists() or not (dest/'INCOMPLETE').exists():
        raise ValueError('unexpected destination state')
    source = derived/'GRCh38_acgtn.fa.gz'
    expected = json.loads((derived/'transformation.json').read_text())['summary']['derived_gz_sha256']
    plain, hashes = verify(source, scratch, dest, expected)
    builder_hash = sha256(BUILDER)
    if builder_hash != '5987ad90e0c4a061cbf4aac7d6b087a284ccf7fe0f9e1835219c567c5e8d1ad2':
        raise ValueError('builder pin mismatch')
    manifest = dict(format='sbds-rebuilt-index-completion-v1', assembly='GRCh38',
        status='WRITTEN_AWAITING_CROSS_NODE_RECHECK', destination=str(dest),
        build_host=socket.gethostname(), source_fasta=str(source),
        source_compressed_sha256=expected, source_plain_sha256=plain,
        builder_sha256=builder_hash, source_override=True, components_sha256=hashes,
        timings={}, build_time=(dest/'build.time').read_text(),
        total_component_bytes=sum((dest/('GRCh38.'+s)).stat().st_size for s in SUFFIXES),
        original_index_untouched=True, recovered_from_job='2875899',
        recovery_job=os.environ['SLURM_JOB_ID'],
        integrity='source and scratch/staged bytes freshly verified; cross-node recheck required')
    with (dest/'completion.json').open('x') as f:
        json.dump(manifest, f, indent=2); f.flush(); os.fsync(f.fileno())
    (dest/'build-host.txt').write_text('BUILD_HOST='+socket.gethostname()+'\n')
    print(json.dumps(manifest))


if __name__ == '__main__':
    main()
