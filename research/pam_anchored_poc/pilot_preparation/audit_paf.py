#!/usr/bin/env python3
"""Independently replay every reported PAF alignment; does not measure recall."""
import argparse, json, re, sys
from pathlib import Path
repo=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(repo/'scripts/release_validation'))
from audit_emitted import independent_cfd, matrices
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--fai',required=True);p.add_argument('--guide',required=True);p.add_argument('--paf',required=True);p.add_argument('--expect-contig');p.add_argument('--expect-strand',choices=('+','-'));p.add_argument('--expect-anchor',type=int);a=p.parse_args()
assert all(x is None for x in (a.expect_contig,a.expect_strand,a.expect_anchor)) or all(x is not None for x in (a.expect_contig,a.expect_strand,a.expect_anchor)), 'expected locus requires contig, strand and anchor'
assert re.fullmatch('[ACGT]{20}',a.guide), 'guide must be exactly 20 A/C/G/T bases'
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
mmatrix,pmatrix=matrices(repo);seen=set();count=scored=0;expected_exact=False
with open(a.paf) as paf:
 for line in paf:
  f=line.rstrip('\n').split('\t');assert len(f)>=12 and f[0]=='pilot' and f[4] in '+-'
  c,s=f[5],f[4];lo,hi=int(f[7]),int(f[8]);assert c in rows and 0<=lo<hi<=rows[c][0]
  assert int(f[1])==len(a.guide) and int(f[2])==0 and int(f[3])==len(a.guide) and int(f[6])==rows[c][0]
  key=(c,s,hi if s=='+' else lo-3);assert key not in seen,'duplicate PAM anchor';seen.add(key)
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
    if op=='I':ga.append(a.guide[qi]);ta.append('-');qi+=1
    elif op=='D':ga.append('-');ta.append(target[ti]);ti+=1
    else:
     same=a.guide[qi]==target[ti];assert same==(op=='=');matches+=same;mismatches+=not same;ga.append(a.guide[qi]);ta.append(target[ti]);qi+=1;ti+=1
  assert qi==len(a.guide) and ti==len(target) and mismatches<=5 and events<=1 and maxgap<=2 and matches/len(a.guide)>=.75
  assert int(f[9])==matches and int(f[10])==sum(n for n,_ in ops) and int(tags['nm'])==mismatches and int(tags['ng'])==events and int(tags['bs'])==maxgap
  if (c,s,key[2])==(a.expect_contig,a.expect_strand,a.expect_anchor):
   assert mismatches==0 and events==0 and cigar=='20=', 'expected locus is not an exact match'
   expected_exact=True
  score=independent_cfd(''.join(ga),''.join(ta),ops,pam,mmatrix,pmatrix);assert ('cf' in tags)==(score is not None)
  if score is not None:assert abs(float(tags['cf'])-score)<.000051;scored+=1
  count+=1
if a.expect_contig is not None:assert expected_exact, 'expected exact target missing'
print(json.dumps({'known_exact_target_recovered':expected_exact if a.expect_contig is not None else None,'status':'pass','reported':count,'duplicate_anchors':0,'independently_scored':scored,'scope':'eligibility of reported sites only; recall not assessed'},sort_keys=True))
