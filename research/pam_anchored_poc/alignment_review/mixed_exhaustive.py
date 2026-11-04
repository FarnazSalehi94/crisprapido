#!/usr/bin/env python3
"""Exhaustive bounded one-event I/D paths against indexed anchor output."""
import csv,functools,itertools,json,subprocess,sys
from pathlib import Path
builder,poc,out=sys.argv[1:];p=Path(out);p.mkdir(exist_ok=False);guide='ACAC';dna='ACGT';rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];refs={}
for n in range(1,6):
 for x in itertools.product('ACN',repeat=n):
  t=''.join(x);name='r'+str(n)+'_'+t;refs[name+'_p']=t+'AGG';refs[name+'_m']=rc(t+'AGG')
(p/'ref.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in refs.items()));(p/'guides.fa').write_text('>g\n'+guide+'\n');subprocess.run([builder,'-r',str(p/'index'),'-f',str(p/'ref.fa')],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
@functools.lru_cache(None)
def eligible(t,m,b,z,frac):
 @functools.lru_cache(None)
 def rec(i,j,mm,matches,state,gaplen):
  if mm>m:return False
  if i==len(guide) and j==len(t):return matches>=1 and matches/len(guide)>=frac
  if i<len(guide) and j<len(t) and rec(i+1,j+1,mm+(guide[i]!=t[j]),matches+(guide[i]==t[j]),2 if state else 0,gaplen):return True
  if b and state!=2 and gaplen<z:
   if i<len(guide) and rec(i+1,j,mm,matches,1,gaplen+1):return True
   if j<len(t) and rec(i,j+1,mm,matches,1,gaplen+1):return True
  return False
 return rec(0,0,0,0,0,0)
results=[]
for m,b,z,f in [(0,0,0,.75),(0,1,2,.75),(1,1,2,.75),(2,1,2,.5),(0,1,1,.5),(4,1,2,0)]:
 truth=set()
 for contig,s in refs.items():
  for strand in '+-':
   ori=s if strand=='+' else rc(s)
   for e in range(1,len(ori)-2):
    if ori[e] not in dna or ori[e+1:e+3]!='GG':continue
    for n in range(max(1,len(guide)-z),len(guide)+z+1):
     if e>=n and eligible(ori[e-n:e],m,b,z,f):truth.add(('g',contig,strand,e if strand=='+' else len(s)-e-3))
 r=subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(p/'guides.fa'),str(m),str(z),str(f),'GG',str(b)],capture_output=True,text=True,check=True);lines=r.stdout.splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'));rows=list(csv.DictReader(lines[i:],delimiter='\t'));got={(x['guide'],x['contig'],x['strand'],int(x['anchor'])-(3 if x['strand']=='-' else 0)) for x in rows}
 result={'m':m,'b':b,'z':z,'f':f,'expected':len(truth),'recovered':len(truth&got),'missed':len(truth-got),'unexpected':len(got-truth),'duplicates':len(rows)-len(got)};results.append(result)
 if got!=truth or len(rows)!=len(got):
  result['missing_examples']=list(sorted(truth-got))[:10];result['unexpected_examples']=list(sorted(got-truth))[:10];(p/'summary.json').write_text(json.dumps(results,indent=2)+'\n');raise AssertionError(result)
(p/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps({'contigs':len(refs),'configurations':len(results),'results':results},sort_keys=True))
