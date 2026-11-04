#!/usr/bin/env python3
"""Deterministic small-reference differential test against direct gap enumeration."""
import collections,csv,json,random,subprocess,sys
from pathlib import Path
R=random.Random(44751); dna='ACGT';rc=lambda x:x.translate(str.maketrans('ACGT','TGCA'))[::-1]
builder,poc,out=sys.argv[1:];p=Path(out);p.mkdir(exist_ok=False); refs={};guides={}
for i in range(16):
    s=''.join(R.choice(dna) for _ in range(65));ori=s if i%2==0 else rc(s);ori=ori[:35]+'AGG'+ori[38:];s=ori if i%2==0 else rc(ori);refs['r'+str(i)]=s
    g=list(ori[15:35]);
    for j in R.sample(range(20),i%6):g[j]=R.choice(dna.replace(g[j],''))
    guides['g'+str(i)]=''.join(g)
(p/'ref.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in refs.items()));(p/'guides.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in guides.items()))
subprocess.run([builder,'-r',str(p/'index'),'-f',str(p/'ref.fa')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
with (p/'out.tsv').open('w') as f,(p/'out.log').open('w') as e:subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(p/'guides.fa'),'5','2','.75'],check=True,stdout=f,stderr=e)
truth=set()
for c,seq in refs.items():
 for strand in '+-':
  ori=seq if strand=='+' else rc(seq)
  for e in range(18,len(seq)-2):
   if ori[e+1:e+3]!='GG':continue
   for delta in (-2,-1,0,1,2):
    n=20+delta
    if e<n:continue
    t=ori[e-n:e]
    for gid,g in guides.items():
     if delta==0:variants=[(g,t)]
     elif delta<0:variants=((g[:i]+g[i-delta:],t) for i in range(n+1))
     else:variants=((g,t[:i]+t[i+delta:]) for i in range(21))
     if any(sum(a!=b for a,b in zip(q,r))<=min(5,len(q)-15) for q,r in variants):truth.add((gid,c,strand,e if strand=='+' else len(seq)-e-3))
ls=(p/'out.tsv').read_text().splitlines();i=next(i for i,x in enumerate(ls) if x.startswith('guide\tcontig\t'));rows=list(csv.DictReader(ls[i:],delimiter='\t'));got={(r['guide'],r['contig'],r['strand'],int(r['anchor'])-(3 if r['strand']=='-' else 0)) for r in rows}
print(json.dumps({'seed':44751,'contigs':len(refs),'guides':len(guides),'expected':len(truth),'recovered':len(truth&got),'missed':len(truth-got),'unexpected':len(got-truth),'duplicates':len(rows)-len(got)},sort_keys=True));assert got==truth and len(got)==len(rows)
