#!/usr/bin/env python3
"""Check whether any missed N-panel anchor has an eligible resolved target."""
import json,sys,collections
from pathlib import Path
from index_check import fasta
p=Path(sys.argv[1]); refs=dict(fasta(p/'ref.fa')); guides=dict(fasta(p/'guides.fa')); rc=lambda x:x.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
missing=json.loads((p/'comparison.json').read_text())['prototype']['missing_keys']; cases=[]
for gid,c,s,a in missing:
 seq=refs[c];g=guides[gid]; ori=seq if s=='+' else rc(seq); end=a if s=='+' else len(seq)-a-3
 resolved=ambiguous=0
 for delta in (-2,-1,0,1,2):
  length=20+delta
  if end<length:continue
  t=ori[end-length:end]
  if delta==0:v=[(g,t)]
  elif delta<0:v=((g[:i]+g[i-delta:],t) for i in range(length+1))
  else:v=((g,t[:i]+t[i+delta:]) for i in range(21))
  for q,r in v:
   mm=sum(x!=y for x,y in zip(q,r))
   if mm<=5 and len(q)-mm>=15:
    if 'N' in t: ambiguous+=1
    else:resolved+=1
 cases.append({'key':[gid,c,s,a],'eligible_resolved_alignments':resolved,'eligible_ambiguous_alignments':ambiguous})
print(json.dumps({'missed':len(cases),'all_have_ambiguous_alignment':all(x['eligible_ambiguous_alignments'] for x in cases),'with_resolved_alternative':sum(x['eligible_resolved_alignments']>0 for x in cases),'cases':cases},indent=2))
