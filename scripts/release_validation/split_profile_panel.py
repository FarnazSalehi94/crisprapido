#!/usr/bin/env python3
"""Split held-out synthetic panel into repetitive and nonrepetitive fixtures."""
import argparse
from pathlib import Path
from resolve_bulge_discrepancies import fasta
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args()
src=Path(a.input);out=Path(a.output);out.mkdir(exist_ok=False)
refs=fasta(src/'ref.fa');guides=fasta(src/'guides.fa')
for label in ('repetitive','nonrepetitive'):
 d=out/label;d.mkdir();isrep=label=='repetitive'
 selected={k:v for k,v in refs.items() if k.startswith('held_repeat_')==isrep}
 queries={k:v for k,v in guides.items() if (k=='held_repeat')==isrep}
 (d/'ref.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in selected.items()))
 (d/'guides.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in queries.items()))
