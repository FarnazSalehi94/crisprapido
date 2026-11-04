#!/usr/bin/env python3
"""Independently replay integrated PAF and enumerate bounded PAM-anchored sites."""
import argparse,json,re,sys
from pathlib import Path
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'scripts/release_validation'));from audit_emitted import independent_cfd,matrices
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--guides',required=True);p.add_argument('--baseline',required=True);p.add_argument('--experimental',required=True);p.add_argument('--out',required=True);p.add_argument('--oracle',action='store_true');a=p.parse_args()
def fasta(path):
 d={};name=None
 for line in Path(path).read_text().splitlines():
  if line.startswith('>'):name=line[1:].split()[0];assert name not in d;d[name]=''
  else:d[name]+=''.join(line.split()).upper()
 return d
refs=fasta(a.reference);guides=fasta(a.guides);mm,pm=matrices(root);rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def replay(path):
 out={};dups=[]
 for line in Path(path).read_text().splitlines():
  f=line.split('\t');assert len(f)>=12
  gid,c,s=f[0],f[5],f[4];g=guides[gid];seq=refs[c];lo,hi=int(f[7]),int(f[8]);assert s in '+-' and 0<=lo<hi<=len(seq)
  target=seq[lo:hi] if s=='+' else rc(seq[lo:hi]);pam=seq[hi:hi+3] if s=='+' else rc(seq[lo-3:lo]);assert len(pam)==3 and pam[0] in 'ACGT' and pam[1:]=='GG'
  tags={x.split(':',2)[0]:x.split(':',2)[2] for x in f[12:]};cigar=tags['cg'];ops=[(int(n),op) for n,op in re.findall(r'(\d+)([=XID])',cigar)];assert ''.join(str(n)+op for n,op in ops)==cigar
  qi=ti=matches=mis=events=gaplen=maxgap=0;prev='';alignedg=[];alignedt=[]
  for n,op in ops:
   if op in 'ID':
    if prev not in ('I','D'):events+=1;gaplen=0
    gaplen+=n;maxgap=max(maxgap,gaplen)
   else:gaplen=0
   prev=op
   for _ in range(n):
    if op=='I':alignedg.append(g[qi]);alignedt.append('-');qi+=1
    elif op=='D':alignedg.append('-');alignedt.append(target[ti]);ti+=1
    else:
     same=g[qi]==target[ti];assert same==(op=='=');matches+=same;mis+=not same;alignedg.append(g[qi]);alignedt.append(target[ti]);qi+=1;ti+=1
  assert qi==len(g) and ti==len(target) and mis<=5 and events<=1 and maxgap<=2 and matches/len(g)>=.75
  assert int(f[1])==len(g) and int(f[2])==0 and int(f[3])==len(g) and int(f[6])==len(seq) and int(f[9])==matches and int(f[10])==sum(n for n,_ in ops)
  assert int(tags['nm'])==mis and int(tags['ng'])==events and int(tags['bs'])==maxgap
  score=independent_cfd(''.join(alignedg),''.join(alignedt),ops,pam,mm,pm);assert ('cf' in tags)==(score is not None)
  if score is not None:assert abs(float(tags['cf'])-score)<.000051
  key=(gid,c,s,hi if s=='+' else lo);row={'start':lo,'end':hi,'cigar':cigar,'mismatches':mis,'bulge_events':events,'max_bulge':maxgap,'matches':matches,'pam':pam,'cfd':score}
  if key in out:dups.append(key)
  out[key]=row
 return out,dups
base,bd=replay(a.baseline);exp,ed=replay(a.experimental);common=set(base)&set(exp);diff=[{'key':k,'baseline':base[k],'experimental':exp[k]} for k in sorted(common) if (base[k]['start'],base[k]['end'],base[k]['cigar'])!=(exp[k]['start'],exp[k]['end'],exp[k]['cigar'])]
def different(x,y):return (x is None)!=(y is None) or (x is not None and y is not None and abs(x-y)>1e-9)
summary={'baseline':len(base),'experimental':len(exp),'both':len(common),'baseline_only':[{'key':k,'row':base[k]} for k in sorted(set(base)-set(exp))],'experimental_only':[{'key':k,'row':exp[k]} for k in sorted(set(exp)-set(base))],'baseline_duplicates':len(bd),'experimental_duplicates':len(ed),'representative_differences':len(diff),'cfd_differences':sum(different(base[k]['cfd'],exp[k]['cfd']) for k in common),'alignments':diff}
def eligible(ori,end,g):
 for ins,dele in [(0,0),(1,0),(2,0),(0,1),(0,2),(1,1)]:
  span=len(g)-ins+dele
  if span<=0 or end<span: continue
  t=ori[end-span:end]
  for k in range(len(g)-ins+1):
   m=sum(x!=y for x,y in zip(g[:k],t[:k]))+sum(x!=y for x,y in zip(g[k+ins:],t[k+dele:]))
   matches=len(g)-ins-m
   if m<=5 and matches>=.75*len(g):return {'guide_gap_bases':ins,'reference_gap_bases':dele,'gap_column':k,'mismatches':m,'matches':matches}
 return None
summary['disagreement_oracle'] = []
for key in sorted(set(base)^set(exp)):
 gid,c,strand,anchor=key;seq=refs[c];ori=seq if strand=='+' else rc(seq);end=anchor if strand=='+' else len(seq)-anchor
 valid=ori[end:end+1] in ('A','C','G','T') and ori[end+1:end+3]=='GG'
 summary['disagreement_oracle'].append({'key':key,'pam_valid':valid,'eligible_witness':eligible(ori,end,guides[gid]) if valid else None})
if a.oracle:
 # All gap events of total length <= 2, including adjacent I/D. At each PAM,
 # independently enumerate the aligned prefix and suffix around the gap.
 truth=set();events=[(0,0),(1,0),(2,0),(0,1),(0,2),(1,1)]
 for gid,g in guides.items():
  for c,seq in refs.items():
   for strand in '+-':
    ori=seq if strand=='+' else rc(seq)
    for end in range(1,len(ori)-2):
     if ori[end] not in 'ACGT' or ori[end+1:end+3]!='GG':continue
     found=False
     for ins,dele in events:
      span=len(g)-ins+dele
      if span<=0 or end<span:continue
      t=ori[end-span:end]
      for k in range(len(g)-ins+1):
       left=g[:k];right=g[k+ins:];m=sum(x!=y for x,y in zip(left,t[:k]))+sum(x!=y for x,y in zip(right,t[k+dele:]));matches=len(g)-ins-m
       if m<=5 and matches>=.75*len(g):found=True;break
      if found:break
     if found:truth.add((gid,c,strand,end if strand=='+' else len(seq)-end))
 summary['oracle']={'expected':len(truth),'baseline_recovered':len(truth&set(base)),'baseline_missing':[k for k in sorted(truth-set(base))],'baseline_unexpected':[k for k in sorted(set(base)-truth)],'experimental_recovered':len(truth&set(exp)),'experimental_missing':[k for k in sorted(truth-set(exp))],'experimental_unexpected':[k for k in sorted(set(exp)-truth)]}
Path(a.out).write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='alignments'}))
