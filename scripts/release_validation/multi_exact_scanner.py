#!/usr/bin/env python3
"""Stream FASTA and report every strict NGG exact hit for multiple guides.

Coordinates are zero-based half-open; both strands and overlaps are retained.
Memory is proportional to the longest FASTA line and guide set, not genome size.
"""
import argparse
import gzip
import json
from collections import defaultdict

from coverage_oracle import DNA, rc


def read_guides(path):
    guides = {}
    name = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                if not name or name in guides:
                    raise ValueError('empty or duplicate guide ID')
                guides[name] = ''
            elif line:
                if name is None:
                    raise ValueError('sequence before FASTA header')
                guides[name] += line.upper()
    if not guides or any(not s or set(s) - DNA for s in guides.values()):
        raise ValueError('guides must be nonempty A/C/G/T sequences')
    return guides


def scan_fasta(handle, guides):
    patterns = defaultdict(list)
    for name, guide in guides.items():
        if not guide or set(guide) - DNA:
            raise ValueError('guides must be nonempty uppercase A/C/G/T')
        patterns[guide].append((name, '+'))
        patterns[rc(guide)].append((name, '-'))
    lengths = sorted({len(p) for p in patterns})
    keep = max(lengths) + 2
    contig, tail, consumed = None, '', 0
    for line in handle:
        line = line.strip()
        if line.startswith('>'):
            contig, tail, consumed = line[1:].split()[0], '', 0
        elif line:
            if contig is None:
                raise ValueError('sequence before FASTA header')
            chunk = line.upper()
            seq = tail + chunk
            base = consumed - len(tail)
            for length in lengths:
                for i in range(max(0, len(tail) - length - 2), len(seq) - length + 1):
                    for name, strand in patterns.get(seq[i:i+length], ()):
                        start = base + i
                        if strand == '+':
                            if i + length + 3 > len(seq) or i + length + 3 <= len(tail):
                                continue
                            pam = seq[i+length:i+length+3]
                        else:
                            if i < 3 or i + length <= len(tail):
                                continue
                            pam = rc(seq[i-3:i])
                        if set(pam) <= DNA and pam[1:] == 'GG':
                            yield dict(guide=name, contig=contig, start=start,
                                       end=start+length, strand=strand, pam=pam)
            consumed += len(chunk)
            tail = seq[-keep:]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta', required=True)
    p.add_argument('--guides', required=True)
    args = p.parse_args()
    guides = read_guides(args.guides)
    opener = gzip.open if args.fasta.endswith('.gz') else open
    with opener(args.fasta, 'rt') as handle:
        for hit in scan_fasta(handle, guides):
            print(json.dumps(hit, sort_keys=True))

if __name__ == '__main__':
    main()
