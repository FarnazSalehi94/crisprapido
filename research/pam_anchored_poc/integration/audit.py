#!/usr/bin/env python3
"""Audit integrated PAF against independent direct-reference oracle and CFD matrix."""
import json, re, sys
from pathlib import Path
repo=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(repo/'scripts/release_validation'))
from audit_emitted import independent_cfd, matrices
sys.path.insert(0,str(repo/'research/pam_anchored_poc/validation'))
from index_check import fasta
ref, guides, oracle, paf = map(Path, sys.argv[1:])
refs=dict(fasta(ref)); gs=dict(fasta(guides)); mmatrix, pmatrix=matrices(repo)
truth={(r['guide'],r['contig'],r['strand'],r['anchor']) for r in map(json.loads,oracle.open())}
rc=lambda x:x.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
seen=set(); count=0; scored=0
for line in paf.read_text().splitlines():
 f=line.split('\t'); assert len(f)>=12
 g,c,s=f[0],f[5],f[4]; lo,hi=int(f[7]),int(f[8]); seq=refs[c]; guide=gs[g]
 assert 0<=lo<hi<=len(seq) and int(f[1])==len(guide) and int(f[2])==0 and int(f[3])==len(guide) and int(f[6])==len(seq)
 key=(g,c,s,hi if s=='+' else lo-3); assert key not in seen;seen.add(key)
 target=seq[lo:hi] if s=='+' else rc(seq[lo:hi]); pam=seq[hi:hi+3] if s=='+' else rc(seq[lo-3:lo]);assert len(pam)==3 and pam[0] in 'ACGT' and pam[1:]=='GG'
 tags={v.split(':',2)[0]:v.split(':',2)[2] for v in f[12:]}; cigar=tags['cg']; ops=[(int(n),op) for n,op in re.findall(r'(\d+)([=XID])',cigar)];assert ''.join(str(n)+op for n,op in ops)==cigar
 qi=ti=matches=mismatches=events=gaplen=maxgap=0; prev=''; a=[];b=[]
 for n,op in ops:
  if op in 'ID':
   if prev not in ('I','D'):events+=1;gaplen=0
   gaplen+=n;maxgap=max(maxgap,gaplen)
  else: gaplen=0
  prev=op
  for _ in range(n):
   if op=='I': a.append(guide[qi]);b.append('-');qi+=1
   elif op=='D': a.append('-');b.append(target[ti]);ti+=1
   else:
    same=guide[qi]==target[ti];assert same==(op=='=');matches+=same;mismatches+=not same;a.append(guide[qi]);b.append(target[ti]);qi+=1;ti+=1
 assert qi==len(guide) and ti==len(target) and mismatches<=5 and events<=1 and maxgap<=2 and matches/len(guide)>=.75
 assert int(f[9])==matches and int(f[10])==sum(n for n,_ in ops) and int(tags['nm'])==mismatches and int(tags['ng'])==events and int(tags['bs'])==maxgap
 score=independent_cfd(''.join(a),''.join(b),ops,pam,mmatrix,pmatrix)
 assert ('cf' in tags)==(score is not None)
 if score is not None: assert abs(float(tags['cf'])-score)<.000051;scored+=1
 count+=1
result={'expected':len(truth),'reported':count,'recovered':len(truth&seen),'missing':len(truth-seen),'unexpected':len(seen-truth),'duplicates':count-len(seen),'independently_scored':scored}
print(json.dumps(result,sort_keys=True))
assert seen==truth
