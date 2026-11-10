#!/usr/bin/env python3
"""Prepare/reuse a trusted, checksum-pinned FASTA offset index, never skip hashes.

Preparation validates every FASTA line with the original reader. Reuse requires
both the trusted index digest and the trusted reference digest. This amortizes
Python wrapping/alphabet validation, not reference identity checking.
"""
import argparse
import json
from pathlib import Path

from adapter import IndexedFasta, digest


class PinnedIndexedFasta(IndexedFasta):
    def __init__(self, path, index, index_sha256, reference_sha256):
        if not index_sha256 or not reference_sha256:
            raise ValueError('pinned index requires index AND reference SHA256')
        self.path = Path(path)
        if digest(index) != index_sha256:
            raise ValueError('reference index checksum mismatch')
        spec = json.loads(Path(index).read_text())
        if spec['format'] != 'crisprapido-fasta-offsets-v1' or spec['reference_sha256'] != reference_sha256:
            raise ValueError('reference index binding mismatch')
        self.sha256 = digest(path)
        if self.sha256 != reference_sha256 or self.path.stat().st_size != spec['reference_bytes']:
            raise ValueError('reference checksum/size mismatch')
        self.rows = {}
        for name, row in spec['rows'].items():
            if not name or len(row) != 4 or any(type(x) is not int for x in row):
                raise ValueError('invalid reference index row')
            length, offset, bases, width = row
            if not 0 < bases <= length or not bases <= width <= bases+2 or not 0 <= offset < spec['reference_bytes']:
                raise ValueError('invalid reference index geometry')
            self.rows[name] = tuple(row)
        if not self.rows or digest(index) != index_sha256:
            raise ValueError('empty or changed reference index')
        self.file = self.path.open('rb')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path, help='new trusted index; never overwrite')
    a = p.parse_args()
    ref = IndexedFasta(a.reference)
    try:
        if digest(a.reference) != ref.sha256:
            raise ValueError('reference changed during preparation')
        spec = dict(format='crisprapido-fasta-offsets-v1', reference_sha256=ref.sha256,
                    reference_bytes=a.reference.stat().st_size, rows=ref.rows)
        with a.output.open('x') as f:
            json.dump(spec, f, sort_keys=True); f.write('\n')
        print(json.dumps(dict(reference_sha256=ref.sha256, index_sha256=digest(a.output))))
    finally:
        ref.close()


if __name__ == '__main__':
    main()
