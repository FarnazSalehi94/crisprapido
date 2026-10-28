#!/usr/bin/env python3
"""Read-only frozen-input/index preflight. Exclusive output, fail closed, no searches.

Reconstructs tag-21 VANILLA indexed text from each full FASTA with the pinned
Columba normalization. SHA256-pins all compressed FASTAs and index components.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys

SUFFIXES = 'bwt brt rev.brt sa.4 sa.bv.4 txt.bin cct pos meta fsid headerSN.bin sna'.split()


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def metadata(row):
    fa, prefix = Path(row['fasta']), row['index_prefix']
    files = [fa, Path(str(fa)+'.fai'), Path(str(fa)+'.gzi')]
    files += [Path(prefix+'.'+s) for s in SUFFIXES]
    for path in files:
        with path.open('rb') as f:
            if not f.read(1):
                raise ValueError('empty input: '+str(path))
    st = fa.stat()
    if st.st_size != int(row['fasta_bytes']) or st.st_mtime_ns != int(row['fasta_mtime_ns']):
        raise ValueError('FASTA metadata changed: '+str(fa))
    for suffix, key in [('.fai', 'fai_sha256'), ('.gzi', 'gzi_sha256')]:
        if sha256(str(fa)+suffix) != row[key]:
            raise ValueError('frozen index hash changed: '+str(fa)+suffix)
    fai = [(x[0], int(x[1])) for x in (s.split('\t') for s in Path(str(fa)+'.fai').read_text().splitlines())]
    header = []
    for line in Path(prefix+'.headerSN.bin').read_text().splitlines():
        fields = dict(x.split(':', 1) for x in line.split('\t')[1:])
        header.append((fields['SN'], int(fields['LN'])))
    if header != fai:
        raise ValueError('contig dictionary/order mismatch: '+prefix)
    data = Path(prefix+'.pos').read_bytes()
    if len(data) % 4:
        raise ValueError('invalid position array')
    positions = struct.unpack('<'+'I'*(len(data)//4), data)
    if (len(positions) != len(fai)+1 or positions[0] != 0 or
            [b-a for a, b in zip(positions, positions[1:])] != [n for _, n in fai]):
        raise ValueError('contig offset/length mismatch: '+prefix)
    data = Path(prefix+'.sna').read_bytes()
    names, pos = [], 0
    while pos < len(data):
        length = struct.unpack_from('<Q', data, pos)[0]
        pos += 8
        names.append(data[pos:pos+length].decode('ascii'))
        pos += length
    if pos != len(data) or names != [n for n, _ in fai]:
        raise ValueError('sequence-name array mismatch: '+prefix)
    if Path(prefix+'.meta').read_text() != '21\n4\nVANILLA\n':
        raise ValueError('unsupported index format: '+prefix)
    if Path(prefix+'.fsid').read_bytes() != struct.pack('<I', 0):
        raise ValueError('expected a single source FASTA: '+prefix)
    return dict(assembly=row['assembly_id'], status='METADATA_PASS', contigs=len(fai),
                bases=sum(n for _, n in fai), files=[str(f) for f in files])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', required=True)
    p.add_argument('--compare-bin', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    out = Path(a.output)
    out.mkdir(exist_ok=False)
    status = dict(status='RUNNING', assemblies=[])
    try:
        plan = Path(a.plan)
        manifest = json.loads((plan/'release_manifest.json').read_text())
        for name, expected in manifest['inputs'].items():
            if sha256(plan/name) != expected:
                raise ValueError('frozen plan hash mismatch: '+name)
        for value in manifest['artifacts'].values():
            if sha256(value['path']) != value['sha256']:
                raise ValueError('artifact hash mismatch: '+value['path'])
        for path, expected in manifest['external_tools'].items():
            if sha256(path) != expected:
                raise ValueError('external tool hash mismatch: '+path)
        guides = list(csv.DictReader((plan/'guides_100.tsv').open(), delimiter='\t'))
        if (len(guides) != 100 or len({g['guide_id'] for g in guides}) != 100 or
                len({g['sequence'] for g in guides}) != 100 or
                any(len(g['sequence']) != 20 or set(g['sequence'])-set('ACGT') for g in guides)):
            raise ValueError('invalid frozen guide set')
        fasta = (plan/'guides_100.fa').read_text()
        if fasta != ''.join('>{}\n{}\n'.format(g['guide_id'], g['sequence']) for g in guides):
            raise ValueError('guide table/FASTA disagree')
        rows = list(csv.DictReader((plan/'assemblies_12.tsv').open(), delimiter='\t'))
        if len(rows) != 12 or len({r['assembly_id'] for r in rows}) != 12:
            raise ValueError('expected twelve frozen assemblies')
        status['metadata'] = [metadata(row) for row in rows]
        (out/'metadata.json').write_text(json.dumps(status['metadata'], indent=2)+'\n')
        status['compare_binary_sha256'] = sha256(a.compare_bin)
        for row in rows:
            asm = row['assembly_id']
            command = [a.compare_bin, row['fasta'], row['index_prefix']+'.txt.bin']
            print(asm+' full sequence comparison started', flush=True)
            result = subprocess.run(command, text=True, capture_output=True)
            (out/(asm+'.text.stdout')).write_text(result.stdout)
            (out/(asm+'.text.stderr')).write_text(result.stderr)
            if result.returncode:
                raise ValueError(asm+' indexed-text verification failed (read/format/content): '+result.stderr.strip())
            evidence = json.loads(result.stdout)
            if evidence['bases'] != int(row['total_bases']) or evidence['contigs'] != int(row['contigs']):
                raise ValueError(asm+' full sequence count mismatch')
            print(asm+' sequence PASS; hashing source and index components', flush=True)
            paths = [Path(row['fasta']), Path(row['fasta']+'.fai'), Path(row['fasta']+'.gzi')]
            paths += [Path(row['index_prefix']+'.'+s) for s in SUFFIXES]
            hashes = {str(path): sha256(path) for path in paths}
            evidence.update(assembly=asm, command=command, hashes=hashes,
                            build_log=Path(row['index_prefix']).parent.joinpath('index_build.time').read_text())
            (out/(asm+'.provenance.json')).write_text(json.dumps(evidence, indent=2)+'\n')
            status['assemblies'].append(dict(assembly=asm, status='PASS'))
            print(asm+' provenance PASS', flush=True)
        status['status'] = 'PASS'
        status['limitations'] = ['Full transformed text/dictionary identity and component hashes reconstruct source binding; not an exhaustive proof of BWT/suffix-array correctness. Functional search coverage must still be audited.']
        (out/'PASS').write_text('All frozen inputs and twelve source/index text bindings passed.\n')
    except Exception as e:
        status.update(status='FAIL', error=str(e))
        print('PREREQUISITE FAILED: '+str(e), file=sys.stderr, flush=True)
    (out/'summary.json').write_text(json.dumps(status, indent=2)+'\n')
    return 0 if status['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
