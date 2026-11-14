#!/usr/bin/env python3
"""Independent post-promotion gate; no search, rebuild, or deletion.

Checks real cross-host evidence, full live component hashes, immutable binding,
source compressed/plain identities, and independent normalized indexed-text
comparison. Produces an exclusive validation record. Never edits input files.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import socket
import subprocess
import time

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def require(ok, message):
    if not ok: raise ValueError(message)


def validate(root, assembly, comparator, comparator_sha256, output):
    started = time.monotonic()
    index_dir = root/'indices'/assembly
    completion = index_dir/'completion.json'
    manifest = json.loads(completion.read_text())
    require(manifest['assembly'] == assembly, 'assembly mismatch')
    require(manifest['status'] == 'TRUSTED_CROSS_NODE', 'not trusted')
    cross = manifest['cross_node_recheck']
    require(cross['result'] == 'PASS' and cross['all_components_match'], 'cross check failed')
    build_host = (index_dir/'build-host.txt').read_text().strip().split('=', 1)[1]
    require(bool(build_host) and cross['build_host'] == build_host and cross['host'] != build_host,
            'missing/invalid independent cross-host verification')
    require(not (index_dir/'INCOMPLETE').exists(), 'incomplete index')
    hashes = {s: digest(index_dir/(assembly+'.'+s)) for s in SUFFIXES}
    require(hashes == manifest['components_sha256'], 'live index component drift')
    source = Path(manifest['source_fasta'])
    require(digest(source) == manifest['source_compressed_sha256'], 'source compressed identity drift')
    h = hashlib.sha256()
    with gzip.open(source, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    require(h.hexdigest() == manifest['source_plain_sha256'], 'decompressed source mismatch')
    binding_dir = root/'runs/sbds_ngg_nrr_v1/00_manifest/bindings'/assembly
    binding_path = binding_dir/'index.crispr-binding.json'
    binding = json.loads(binding_path.read_text())
    require(binding['components_sha256'] == hashes, 'binding component identity mismatch')
    require(binding['reference_sha256'] == h.hexdigest(), 'binding/reference mismatch')
    require(digest(binding_dir/'reference.fa') == h.hexdigest(), 'bound reference changed')
    require(not (binding_dir/'INCOMPLETE').exists() and (binding_dir/'binding.json').exists(),
            'binding not committed')
    for s in SUFFIXES:
        require((binding_dir/('index.'+s)).resolve() == (index_dir/(assembly+'.'+s)).resolve(),
                'binding points outside authoritative index')
    require(digest(comparator) == comparator_sha256, 'comparator pin mismatch')
    result = subprocess.run([str(comparator), str(source), str(index_dir/(assembly+'.txt.bin'))],
                            text=True, capture_output=True, timeout=900, check=True)
    comparison = json.loads(result.stdout)
    require(comparison['status'] == 'PASS', 'indexed text comparison failed')
    record = dict(status='PASS', assembly=assembly, host=socket.gethostname(),
                  seconds=time.monotonic()-started, completion_sha256=digest(completion),
                  binding_sha256=digest(binding_path), components_sha256=hashes,
                  independent_text_comparison=comparison,
                  comparator_sha256=comparator_sha256, comparator=str(comparator),
                  source_sha256=manifest['source_compressed_sha256'])
    with output.open('x') as f: json.dump(record, f, indent=2); f.write('\n')
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--assembly', required=True)
    p.add_argument('--comparator', type=Path, required=True)
    p.add_argument('--comparator-sha256', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(validate(a.root, a.assembly, a.comparator, a.comparator_sha256, a.output)))
