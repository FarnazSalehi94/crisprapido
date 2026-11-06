#!/usr/bin/env python3
"""Count Ns in bounded first/last contig windows using a SHA-pinned FAI."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--fai',required=True);p.add_argument('--window',type=int,default=1000);a=p.parse_args();stage=Path(a.stage);fai=Path(a.fai)
assert hashlib.sha256(fai.read_bytes()).hexdigest()==json.loads((stage/'import.json').read_text())['components']['fai']['sha256']
rows=[x.split('\t') for x in fai.read_text().splitlines()];result=[]
with (stage/'reference.fa').open('rb') as f:
 for row in rows:
  name,L,off,b,w=row[:5];L,off,b,w=map(int,(L,off,b,w))
  for side,lo,hi in [('first',0,min(a.window,L)),('last',max(0,L-a.window),L)]:
   seq=[]
   while lo<hi:
    n=min(hi-lo,b-lo%b);f.seek(off+(lo//b)*w+lo%b);chunk=f.read(n);assert len(chunk)==n and set(chunk)<=set(b'ACGTNacgtn');seq.append(chunk.upper());lo+=n
   n=b''.join(seq).count(b'N')
   if n:result.append({'contig':name,'side':side,'count':n})
print(json.dumps({'window':a.window,'contigs':len(rows),'windows_with_N':result,'qualification':'bounded FASTA inspection; no indexed guide search'}))
