#!/usr/bin/env python3
"""Exercise every 2-base PAM and supported one-gap options against direct oracle."""
import csv,itertools,json,re,subprocess,sys
from pathlib import Path
DNA='ACGT';rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
builder,poc,out=sys.argv[1:];p=Path(out);p.mkdir(exist_ok=False);g='ACGTTGCACTGAGTACCGTA';refs={}
for pam in map(''.join,itertools.product(DNA,repeat=2)):
 for n in DNA:
  for strand in '+-':
   # Some targets have mismatch, 1I, 2D, and reference N.
   for mode in range(4):
    t=[g,g[:5]+'T'+g[6:],g[:8]+g[9:],g[:8]+'NN'+g[8:]][mode];s=t+n+pam
    refs['%s_%s_%s_%d'%(pam,n,'p' if strand=='+' else 'm',mode)]=s if strand=='+' else rc(s)
(p/'ref.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in refs.items()));(p/'guides.fa').write_text('>g\n'+g+'\n')
subprocess.run([builder,'-r',str(p/'index'),'-f',str(p/'ref.fa')],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
def oracle(pam,m,b,z,f):
 truth=set()
 for c,seq in refs.items():
  for strand in '+-':
   ori=seq if strand=='+' else rc(seq)
   for e in range(1,len(ori)-2):
    if ori[e] not in DNA or ori[e+1:e+3]!=pam:continue
    for delta in range(-z if b else 0,(z if b else 0)+1):
     n=20+delta
     if n<=0 or e<n:continue
     t=ori[e-n:e]
     if delta==0:
      v=[(g,t)]
      if b and z>=2:v.extend((g[:i]+g[i+1:],t[:i]+t[i+1:]) for i in range(20))
     elif delta<0:v=((g[:i]+g[i-delta:],t) for i in range(n+1))
     else:v=((g,t[:i]+t[i+delta:]) for i in range(21))
     if any((sum(a!=b for a,b in zip(q,r))<=m and len(q)-sum(a!=b for a,b in zip(q,r))>=f*20 and len(q)-sum(a!=b for a,b in zip(q,r))>=1) for q,r in v):truth.add(('g',c,strand,e if strand=='+' else len(seq)-e-3))
 return truth
stats=[]; cfd_rows=[]; cfd_header=None
for pam in map(''.join,itertools.product(DNA,repeat=2)):
 for m,b,z,f in [(0,0,2,.75),(0,1,2,.75),(2,1,1,.75),(5,1,2,.75),(5,1,0,.75),(20,0,0,0.0)]:
  run=subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(p/'guides.fa'),str(m),str(z),str(f),pam,str(b)],text=True,capture_output=True,check=True); lines=run.stdout.splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'));rows=list(csv.DictReader(lines[i:],delimiter='\t')); got={(r['guide'],r['contig'],r['strand'],int(r['anchor'])-(3 if r['strand']=='-' else 0)) for r in rows};want=oracle(pam,m,b,z,f)
  assert got==want and len(got)==len(rows),(pam,m,b,z,len(want-got),len(got-want));
  for row in rows:
   target=refs[row['contig']][int(row['start']):int(row['end'])];target=target if row['strand']=='+' else rc(target);qi=ti=matches=mismatches=ins=dele=events=gaplen=0;previous=''
   ops=re.findall(r'(\d+)([=XID])',row['cigar']);assert ''.join(n+op for n,op in ops)==row['cigar']
   for raw,op in ops:
    n=int(raw)
    if op in 'ID':
     if previous not in 'ID' or previous=='':events+=1;gaplen=0
     gaplen+=n;assert gaplen<=z
    else:gaplen=0
    previous=op
    for _ in range(n):
     if op=='I':qi+=1;ins+=1
     elif op=='D':ti+=1;dele+=1
     else:
      eq=g[qi]==target[ti];assert eq==(op=='=');qi+=1;ti+=1;matches+=eq;mismatches+=not eq
   assert qi==20 and ti==len(target) and mismatches<=m and events<=b and matches>=1 and matches>=20*f
   assert (mismatches,ins,dele,events,matches)==tuple(int(row[k]) for k in ('mismatches','inserted_guide_bases','deleted_guide_bases','bulge_events','match_count'))
   assert abs(float(row['match_fraction'])-matches/20)<1e-6
  if (m,b,z,f)==(5,1,2,.75):
   cfd_header=lines[i];cfd_rows.extend(lines[i+1:])
  stats.append({'pam':pam,'m':m,'b':b,'z':z,'f':f,'anchors':len(got)})
# Unsupported mode and malformed parameters must fail rather than return hits.
for tail in [('5','2','.75','GG','2'),('5','3','.75','GG','1'),('-1','2','.75','GG','1'),('5','2','nan','GG','1'),('5','2','.75','GN','1')]:
 r=subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(p/'guides.fa')]+list(tail),capture_output=True,text=True);assert r.returncode!=0 and not any(l.startswith('g\t') for l in r.stdout.splitlines()),tail
# Invalid guide in any batch position must fail before reporting any rows.
for name,text in [('N','>good\n'+g+'\n>bad\n'+g[:10]+'N'+g[11:]+'\n'),('duplicate','>same\n'+g+'\n>same\n'+g+'\n')]:
 file=p/('bad_'+name+'.fa');file.write_text(text)
 run=subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(file),'5','2','.75','GG','1'],capture_output=True,text=True)
 assert run.returncode!=0 and not any(x.startswith(('good\t','same\t')) for x in run.stdout.splitlines()),name
(p/'cfd_all_pams.tsv').write_text(cfd_header+'\n'+'\n'.join(cfd_rows)+'\n')
(p/'summary.json').write_text(json.dumps({'runs':stats,'unsupported_checked':5},indent=2)+'\n');print(json.dumps({'pam_patterns':16,'option_runs':len(stats),'unsupported_checked':5,'invalid_guide_batches_checked':2,'all_pass':True}))
