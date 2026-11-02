#!/usr/bin/env python3
"""Check every oracle interval for any overlapping e=5/e=0 raw candidate."""
import argparse
import json
from pathlib import Path
import re
from resolve_bulge_discrepancies import fasta,locus,overlaps
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--panel',required=True);p.add_argument('--output',required=True);a=p.parse_args()
d=Path(a.panel);refs=fasta(d/'ref.fa');truth=json.loads((d/'truth.json').read_text());raw=set()
for name in ('edit.sam','exact.sam'):
 for line in (d/name).read_text().splitlines():
  if line.startswith('@'):continue
  f=line.split('\t');flag=int(f[1])
  if flag&4:continue
  start=int(f[3])-1;span=sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]) if op in 'MDN=X')
  raw.add((f[0],f[2],start,start+span,'-' if flag&16 else '+'))
no_overlap=[tuple(x['key']) for x in truth if not any(overlaps(tuple(x['key']),c) for c in raw)]
all_loci={locus(tuple(x['key']),refs) for x in truth}
covered_loci={locus(tuple(x['key']),refs) for x in truth if any(overlaps(tuple(x['key']),c) for c in raw)}
result={'truth_intervals':len(truth),'unique_raw_intervals':len(raw),'intervals_without_overlap':no_overlap,
        'truth_loci':len(all_loci),'loci_without_any_overlap':sorted(all_loci-covered_loci)}
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
