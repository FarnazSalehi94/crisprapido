#!/usr/bin/env python3
"""Check every replacement of three consecutive reference Ns on both strands."""
import argparse,csv,itertools,json,subprocess
from pathlib import Path
DNA='ACGT';rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
a=argparse.ArgumentParser();a.add_argument('--out',required=True);a.add_argument('--builder',required=True);a.add_argument('--poc',required=True);a=a.parse_args();p=Path(a.out);p.mkdir(exist_ok=False)
guides={'complex':'ACGTTGCACTGAGTACCGTA','repeat':'AC'*10};refs={}
for k,g in guides.items():
 for strand in '+-':
  t=g[:7]+'NNN'+g[10:]; s='CC'+t+'AGGGAC';refs[k+('_plus' if strand=='+' else '_minus')]=s if strand=='+' else rc(s)
(p/'ref.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in refs.items()));(p/'guides.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in guides.items()))
truth=set()
for contig,seq in refs.items():
 for strand in '+-':
  ori=seq if strand=='+' else rc(seq)
  for e in range(18,len(ori)-2):
   if ori[e] not in DNA or ori[e+1:e+3]!='GG':continue
   for delta in (-2,-1,0,1,2):
    n=20+delta
    if e<n:continue
    t=ori[e-n:e]
    for gid,g in guides.items():
     if delta==0:v=[(g,t)]
     elif delta<0:v=((g[:i]+g[i-delta:],t) for i in range(n+1))
     else:v=((g,t[:i]+t[i+delta:]) for i in range(21))
     if any(sum(x!=y for x,y in zip(q,r))<=min(5,len(q)-15) for q,r in v):truth.add((gid,contig,strand,e if strand=='+' else len(seq)-e-3))
results=[]
for triplet in itertools.product(DNA,repeat=3):
 txt=''.join(triplet);d=p/txt;d.mkdir();idx=d/'index';f=d/'idx.fa';records={}
 for name,seq in refs.items():
  chars=iter(triplet);records[name]=''.join(next(chars) if c=='N' else c for c in seq)
 f.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in records.items()));subprocess.run([a.builder,'-r',str(idx),'-f',str(f)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
 o=subprocess.run([a.poc,str(idx),str(p/'ref.fa'),str(p/'guides.fa'),'5','2','.75'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=True).stdout.splitlines(); i=next(i for i,x in enumerate(o) if x.startswith('guide\tcontig\t')); rows=list(csv.DictReader(o[i:],delimiter='\t'));found={(x['guide'],x['contig'],x['strand'],int(x['anchor'])-(3 if x['strand']=='-' else 0)) for x in rows}
 results.append({'substitute':txt,'expected':len(truth),'recovered':len(truth&found),'missed':len(truth-found),'unexpected':len(found-truth),'duplicates':len(rows)-len(found)})
 assert found==truth and len(rows)==len(found),results[-1]
summary={'combinations':len(results),'all_pass':all(not x['missed'] and not x['unexpected'] and not x['duplicates'] for x in results),'oracle_anchors':len(truth),'results':results};(p/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'combinations':len(results),'all_pass':summary['all_pass'],'oracle_anchors':len(truth)}))
