#!/usr/bin/env python3
"""Summarize PAM-anchor recovery by oracle category without interval equivalence."""
import argparse
import json
from pathlib import Path
from resolve_bulge_discrepancies import fasta, locus
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--panel',required=True);p.add_argument('--output',required=True);a=p.parse_args()
d=Path(a.panel);refs=fasta(d/'ref.fa');truth=json.loads((d/'truth.json').read_text());reported=set()
for line in (d/'out.paf').read_text().splitlines():
 f=line.split('\t');reported.add(locus((f[0],f[5],int(f[7]),int(f[8]),f[4]),refs))
res=[]
for category in ('exact','mismatch','I','D','mismatch_bulge'):
 expected={locus(tuple(x['key']),refs) for x in truth if category in x['categories']}
 res.append({'category':category,'oracle_loci':len(expected),'recovered':len(expected&reported),
             'recall':len(expected&reported)/len(expected) if expected else None,
             'missed':sorted(expected-reported)})
Path(a.output).write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
