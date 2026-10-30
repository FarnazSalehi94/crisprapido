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


def scan_fasta(handle, guides, chunk_size=1024 * 1024):
    """Scan with C-level substring search, retaining line/contig boundary context."""
    if chunk_size < 1:
        raise ValueError('chunk size must be positive')
    patterns = defaultdict(list)
    for name, guide in guides.items():
        if not guide or set(guide) - DNA:
            raise ValueError('guides must be nonempty uppercase A/C/G/T')
        patterns[guide].append((name, '+'))
        patterns[rc(guide)].append((name, '-'))
    if not patterns:
        raise ValueError('empty guide set')
    keep = max(map(len, patterns)) + 2
    contig, tail, consumed = None, '', 0
    pending, pending_len = [], 0

    def scan_chunk(chunk):
        nonlocal tail, consumed
        seq = tail + chunk
        old = len(tail)
        base = consumed - old
        for pattern, entries in patterns.items():
            start = 0
            while True:
                i = seq.find(pattern, start)
                if i < 0:
                    break
                start = i + 1
                for name, strand in entries:
                    end = i + len(pattern)
                    if strand == '+':
                        if end + 3 > len(seq) or end + 3 <= old:
                            continue
                        pam = seq[end:end+3]
                    else:
                        if i < 3 or end <= old:
                            continue
                        pam = rc(seq[i-3:i])
                    if set(pam) <= DNA and pam[1:] == 'GG':
                        yield dict(guide=name, contig=contig, start=base+i,
                                   end=base+end, strand=strand, pam=pam)
        consumed += len(chunk)
        tail = seq[-keep:]

    for line in handle:
        line = line.strip()
        if line.startswith('>'):
            if pending:
                yield from scan_chunk(''.join(pending))
            contig, tail, consumed = line[1:].split()[0], '', 0
            pending, pending_len = [], 0
        elif line:
            if contig is None:
                raise ValueError('sequence before FASTA header')
            pending.append(line.upper())
            pending_len += len(line)
            if pending_len >= chunk_size:
                yield from scan_chunk(''.join(pending))
                pending, pending_len = [], 0
    if pending:
        yield from scan_chunk(''.join(pending))


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
