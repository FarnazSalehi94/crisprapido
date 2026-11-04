#!/usr/bin/env python3
"""Compare PAM-start keys, independently audit emitted CIGAR and summarize misses."""
import csv,json,re,sys,collections
from pathlib import Path
from index_check import fasta
p=Path(sys.argv[1]); refs=dict(fasta(p/'ref.fa')); guides=dict(fasta(p/'guides.fa')); rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
truth={(x['guide'],x['contig'],x['strand'],x['anchor']):x['categories'] for x in map(json.loads,(p/'oracle.jsonl').open())}
def key(g,c,s,a):return (g,c,s,a-(3 if s=='-' else 0))
def validate(g,c,s,b,e,cigar):
 assert 0<=b<e<=len(refs[c]); t=refs[c][b:e]; pam=refs[c][e:e+3] if s=='+' else rc(refs[c][b-3:b]);
 if s=='-': assert b>=3;t=rc(t)
 assert len(pam)==3 and pam[0] in 'ACGT' and pam[1:]=='GG'
 ops=re.findall(r'(\d+)([=XID])',cigar);assert ''.join(n+o for n,o in ops)==cigar
 qi=ti=mm=match=gaps=ins=dele=0
 for n,op in ops:
  n=int(n)
  if op in 'ID': gaps+=1;assert n<=2
  for _ in range(n):
   if op=='I':qi+=1;ins+=1
   elif op=='D':ti+=1;dele+=1
   else:
    assert qi<len(guides[g]) and ti<len(t); eq=guides[g][qi]==t[ti];assert eq==(op=='=');mm+=not eq;match+=eq;qi+=1;ti+=1
 assert qi==len(guides[g]) and ti==len(t) and mm<=5 and gaps<=1 and match>=15, (g,c,s,b,e,cigar,mm,gaps,match)
 return {'mismatches':mm,'gaps':gaps,'inserted_guide_bases':ins,'deleted_guide_bases':dele,'target_has_N':'N' in t}
lines=(p/'poc.tsv').read_text().splitlines();idx=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t')); rows=list(csv.DictReader(lines[idx:],delimiter='\t'))
proto=[]
for x in rows:
 g=x['guide'];c=x['contig'];s=x['strand'];b=int(x['start']);e=int(x['end']);a=int(x['anchor']);assert a==(e if s=='+' else b)
 v=validate(g,c,s,b,e,x['cigar']); assert all(v[k]==int(x[k]) for k in ('mismatches','inserted_guide_bases','deleted_guide_bases'));proto.append(key(g,c,s,a))
base=[]
for line in (p/'baseline.paf').read_text().splitlines():
 f=line.split('\t');g,c,s=f[0],f[5],f[4];b,e=int(f[7]),int(f[8]); cig=next(x[5:] for x in f[12:] if x.startswith('cg:Z:'));validate(g,c,s,b,e,cig);base.append(key(g,c,s,e if s=='+' else b))
def summarize(x):
 s=set(x);missing=set(truth)-s;extra=s-set(truth)
 return {'reported':len(x),'unique':len(s),'expected':len(truth),'recovered':len(set(truth)&s),'missed':len(missing),'unexpected':len(extra),'duplicates':len(x)-len(s),'missing_keys':[list(k) for k in sorted(missing)],'extra_keys':[list(k) for k in sorted(extra)]}
r={'prototype':summarize(proto),'baseline':summarize(base),'oracle_categories':dict(collections.Counter(c for cats in truth.values() for c in cats)),'index':json.loads((p/'index-check.json').read_text())}
print(json.dumps(r,indent=2,sort_keys=True))
