#!/usr/bin/env python3
"""Audit bounded benchmark alignments without an exhaustive oracle."""
import csv,re,sys
from pathlib import Path
from index_check import fasta
ref,guide,poc,base=sys.argv[1:];refs=dict(fasta(ref));guides=dict(fasta(guide));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def audit(g,c,s,b,e,cigar):
 seq=refs[c];assert 0<=b<e<=len(seq);pam=seq[e:e+3] if s=='+' else rc(seq[b-3:b]); assert len(pam)==3 and pam[0] in 'ACGT' and pam[1:]=='GG'
 t=seq[b:e] if s=='+' else rc(seq[b:e]);q=guides[g];ops=re.findall(r'(\d+)([=XID])',cigar);assert ''.join(n+o for n,o in ops)==cigar
 qi=ti=mm=matches=gaps=0
 for n,op in ops:
  n=int(n)
  if op in 'ID':gaps+=1;assert n<=2
  for _ in range(n):
   if op=='I':qi+=1
   elif op=='D':ti+=1
   else:
    assert qi<len(q) and ti<len(t);eq=q[qi]==t[ti];assert eq==(op=='=');qi+=1;ti+=1;mm+=not eq;matches+=eq
 assert qi==len(q) and ti==len(t) and mm<=5 and gaps<=1 and matches>=.75*len(q)
 return (g,c,s,e if s=='+' else b)
ls=Path(poc).read_text().splitlines();i=next(i for i,x in enumerate(ls) if x.startswith('guide\tcontig\t'));rows=list(csv.DictReader(ls[i:],delimiter='\t')); a=[]
for r in rows:
 b,e=int(r['start']),int(r['end']);key=audit(r['guide'],r['contig'],r['strand'],b,e,r['cigar']);assert int(r['anchor'])==key[-1];a.append(key)
b=[]
for line in Path(base).read_text().splitlines():
 f=line.split('\t');cig=next(x[5:] for x in f[12:] if x.startswith('cg:Z:'));b.append(audit(f[0],f[5],f[4],int(f[7]),int(f[8]),cig))
assert len(a)==len(set(a)) and len(b)==len(set(b)) and set(a)==set(b)
print('audited prototype/baseline anchors:',len(a),len(b),'matching:',len(set(a)&set(b)))
