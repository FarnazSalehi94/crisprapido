#!/usr/bin/env python3
"""Independently replay every reported PAF alignment; does not measure recall."""
import argparse, json, re, sys
from pathlib import Path
repo=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(repo/'scripts/release_validation'))
from audit_emitted import independent_cfd, matrices
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--fai',required=True);p.add_argument('--guides',required=True);p.add_argument('--targets',required=True);p.add_argument('--paf',required=True);p.add_argument('--expect-contig');p.add_argument('--expect-strand',choices=('+','-'));p.add_argument('--expect-anchor',type=int);a=p.parse_args()
guides={};name=None
for line in Path(a.guides).read_text().splitlines():
 if line.startswith('>'):
  name=line[1:].split()[0];assert name and name not in guides;guides[name]=''
 else:
  assert name is not None;guides[name]+=line.strip().upper()
assert guides and all(re.fullmatch('[ACGT]{20}',g) for g in guides.values())
targets=json.loads(Path(a.targets).read_text());assert set(targets)==set(guides)
found=set()
per_guide={q:{'reported':0,'independently_scored':0,'duplicate_anchors':0} for q in guides}
rows={}
for line in open(a.fai):
 f=line.rstrip().split('\t'); assert len(f)>=5 and f[0] not in rows
 rows[f[0]]=tuple(map(int,f[1:5]))
ref=open(a.reference,'rb')
def fetch(name,lo,hi):
 length,off,bases,width=rows[name];assert 0<=lo<=hi<=length and bases>0 and width>=bases
 out=[]
 while lo<hi:
  n=min(hi-lo,bases-lo%bases);ref.seek(off+(lo//bases)*width+lo%bases);x=ref.read(n).decode('ascii').upper();assert len(x)==n and set(x)<=set('ACGTN');out.append(x);lo+=n
 return ''.join(out)
rc=lambda x:x.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
mmatrix,pmatrix=matrices(repo);seen=set();count=scored=0
with open(a.paf) as paf:
 for line in paf:
  f=line.rstrip('\n').split('\t');assert len(f)>=12 and f[0] in guides and f[4] in '+-'
  guide=guides[f[0]]
  c,s=f[5],f[4];lo,hi=int(f[7]),int(f[8]);assert c in rows and 0<=lo<hi<=rows[c][0]
  assert int(f[1])==len(guide) and int(f[2])==0 and int(f[3])==len(guide) and int(f[6])==rows[c][0]
  key=(f[0],c,s,hi if s=='+' else lo-3);assert key not in seen,'duplicate PAM anchor';seen.add(key)
  target=fetch(c,lo,hi);target=target if s=='+' else rc(target)
  pam=fetch(c,hi,hi+3) if s=='+' else rc(fetch(c,lo-3,lo));assert pam[0] in 'ACGT' and pam[1:]=='GG'
  tags={v.split(':',2)[0]:v.split(':',2)[2] for v in f[12:]};cigar=tags['cg'];ops=[(int(n),op) for n,op in re.findall(r'(\d+)([=XID])',cigar)];assert ''.join(str(n)+op for n,op in ops)==cigar
  qi=ti=matches=mismatches=events=gaplen=maxgap=0;prev='';ga=[];ta=[]
  for n,op in ops:
   assert n>0
   if op in 'ID':
    if prev not in ('I','D'):events+=1;gaplen=0
    gaplen+=n;maxgap=max(maxgap,gaplen)
   else:gaplen=0
   prev=op
   for _ in range(n):
    if op=='I':ga.append(guide[qi]);ta.append('-');qi+=1
    elif op=='D':ga.append('-');ta.append(target[ti]);ti+=1
    else:
     same=guide[qi]==target[ti];assert same==(op=='=');matches+=same;mismatches+=not same;ga.append(guide[qi]);ta.append(target[ti]);qi+=1;ti+=1
  assert qi==len(guide) and ti==len(target) and mismatches<=5 and events<=1 and maxgap<=2 and matches/len(guide)>=.75
  assert int(f[9])==matches and int(f[10])==sum(n for n,_ in ops) and int(tags['nm'])==mismatches and int(tags['ng'])==events and int(tags['bs'])==maxgap
  expected=targets[f[0]]
  if (c,s,key[3])==(expected['contig'],expected['strand'],expected['anchor']):
   assert mismatches==0 and events==0 and cigar=='20=', 'expected locus is not an exact match'
   found.add(f[0])
  score=independent_cfd(''.join(ga),''.join(ta),ops,pam,mmatrix,pmatrix);assert ('cf' in tags)==(score is not None)
  if score is not None:assert abs(float(tags['cf'])-score)<.000051;scored+=1;per_guide[f[0]]['independently_scored']+=1
  count+=1;per_guide[f[0]]['reported']+=1
assert found==set(guides), 'expected exact target missing for '+str(set(guides)-found)
print(json.dumps({'known_exact_targets_recovered':sorted(found),'per_guide':per_guide,'status':'pass','reported':count,'duplicate_anchors':0,'independently_scored':scored,'scope':'eligibility of reported sites only; recall not assessed'},sort_keys=True))
