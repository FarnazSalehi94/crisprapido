#!/usr/bin/env python3
"""Compare streaming exact truth with real Columba exact SAM, checking FASTA/PAM."""
import argparse
import json
from pathlib import Path
import subprocess
from multi_exact_scanner import read_guides
from coverage_oracle import rc

p = argparse.ArgumentParser(description=__doc__)
for name in ('oracle', 'sam', 'fasta', 'guides', 'samtools', 'output'):
    p.add_argument('--'+name, required=True)
a = p.parse_args()
guides = read_guides(a.guides)
truth_rows = [json.loads(x) for x in Path(a.oracle).read_text().splitlines()]
def key(x):
    return (x['guide'], x['contig'], x['start'], x['end'], x['strand'])
truth = {key(x) for x in truth_rows}
if len(truth) != len(truth_rows):
    raise ValueError('duplicate oracle exact loci')
found, raw, rejected = set(), 0, []
for line in Path(a.sam).read_text().splitlines():
    if line.startswith('@'): continue
    f = line.split('\t'); flag = int(f[1])
    if flag & 4: continue
    raw += 1
    if f[0] not in guides: raise ValueError('unknown guide')
    guide = guides[f[0]]; start = int(f[3])-1; end = start + len(guide)
    if f[5] != str(len(guide))+'M':
        raise ValueError('unexpected exact SAM CIGAR: '+f[5])
    strand = '-' if flag & 16 else '+'
    lo, hi = (start, end+3) if strand == '+' else (start-3, end)
    if lo < 0:
        rejected.append({'sam':line,'reason':'incomplete PAM'}); continue
    region = '{}:{}-{}'.format(f[2], lo+1, hi)
    output = subprocess.run([a.samtools, 'faidx', a.fasta, region], check=True,
                            stdout=subprocess.PIPE, universal_newlines=True).stdout
    sequence = ''.join(output.splitlines()[1:]).upper()
    oriented = sequence if strand == '+' else rc(sequence)
    if len(oriented) != len(guide)+3 or set(oriented) - set('ACGT') or oriented[:len(guide)] != guide or oriented[-2:] != 'GG':
        rejected.append({'sam':line,'reason':'reference target or strict PAM fails','fetched_oriented':oriented}); continue
    found.add((f[0], f[2], start, end, strand))
result = {'oracle_loci':len(truth), 'raw_mapped_rows':raw, 'strict_sam_loci':len(found),
          'recovered':len(truth & found), 'recall':len(truth & found)/len(truth) if truth else None,
          'precision':len(truth & found)/len(found) if found else None,
          'missed':sorted(truth-found), 'unexpected':sorted(found-truth), 'rejected_sam_rows':rejected,
          'per_guide':{g: {'truth':sum(k[0]==g for k in truth), 'recovered':sum(k[0]==g for k in truth & found)} for g in guides}}
Path(a.output).write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, sort_keys=True))
