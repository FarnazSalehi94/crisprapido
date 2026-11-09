#!/usr/bin/env python3
"""Time independent replay/CFD and byte-output for already-reported PAF; no search."""
import argparse,hashlib,json,re,sys,tempfile,time
from pathlib import Path
repo=Path(__file__).resolve().parents[3];sys.path.insert(0,str(repo/'scripts/release_validation'));from audit_emitted import independent_cfd,matrices
p=argparse.ArgumentParser();p.add_argument('--paf',required=True);p.add_argument('--reference',required=True);p.add_argument('--fai',required=True);p.add_argument('--guide',required=True);p.add_argument('--out',required=True);a=p.parse_args();timings={}
def mark(name,t):timings[name]=time.perf_counter()-t
start=time.perf_counter();rows={}
for line in open(a.fai):
 f=line.rstrip().split('\t');assert len(f)>=5 and f[0] not in rows;rows[f[0]]=tuple(map(int,f[1:5]))
records=[];rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];raw=Path(a.paf).read_bytes();seqf=open(a.reference,'rb')
def fetch(c,lo,hi):
 length,offset,bases,width=rows[c];assert 0<=lo<=hi<=length;v=[]
 while lo<hi:
  n=min(hi-lo,bases-lo%bases);seqf.seek(offset+(lo//bases)*width+lo%bases);chunk=seqf.read(n).decode('ascii').upper();assert len(chunk)==n;v.append(chunk);lo+=n
 return ''.join(v)
seen=set()
for line in raw.decode('ascii').splitlines():
 f=line.split('\t');assert len(f)>=12 and f[0]=='pilot' and f[4] in '+-';c=f[5];strand=f[4];lo,hi=int(f[7]),int(f[8]);assert c in rows and int(f[6])==rows[c][0]
 key=(c,strand,hi if strand=='+' else lo-3);assert key not in seen;seen.add(key)
 target=fetch(c,lo,hi);target=target if strand=='+' else rc(target);pam=fetch(c,hi,hi+3) if strand=='+' else rc(fetch(c,lo-3,lo));assert pam[0] in 'ACGT' and pam[1:]=='GG'
 tags={v.split(':',2)[0]:v.split(':',2)[2] for v in f[12:]};ops=[(int(n),op) for n,op in re.findall(r'(\d+)([=XID])',tags['cg'])];assert ''.join(str(n)+op for n,op in ops)==tags['cg'];records.append((f,target,pam,tags,ops))
mark('parse_and_reference_random_access_seconds',start)
start=time.perf_counter();aligned=[]
for f,target,pam,tags,ops in records:
 qi=ti=matches=mm=events=gap=maxgap=0;prev='';ga=[];ta=[]
 for n,op in ops:
  if op in 'ID':
   if prev not in ('I','D'):events+=1;gap=0
   gap+=n;maxgap=max(maxgap,gap)
  else:gap=0
  prev=op
  for _ in range(n):
   if op=='I':ga.append(a.guide[qi]);ta.append('-');qi+=1
   elif op=='D':ga.append('-');ta.append(target[ti]);ti+=1
   else:
    same=a.guide[qi]==target[ti];assert same==(op=='=');matches+=same;mm+=not same;ga.append(a.guide[qi]);ta.append(target[ti]);qi+=1;ti+=1
 assert qi==len(a.guide) and ti==len(target) and mm<=5 and events<=1 and maxgap<=2 and matches/len(a.guide)>=.75
 assert (matches,mm,events,maxgap)==(int(f[9]),int(tags['nm']),int(tags['ng']),int(tags['bs']))
 aligned.append((''.join(ga),''.join(ta),ops,pam,tags))
mark('independent_emitted_alignment_replay_seconds',start)
start=time.perf_counter();mmatrix,pmatrix=matrices(repo);scored=0
for guide,target,ops,pam,tags in aligned:
 score=independent_cfd(guide,target,ops,pam,mmatrix,pmatrix);assert ('cf' in tags)==(score is not None)
 if score is not None:assert abs(float(tags['cf'])-score)<.000051;scored+=1
mark('independent_cfd_scoring_seconds_including_matrix_load',start)
start=time.perf_counter()
with tempfile.TemporaryFile() as f:
 f.write(raw);f.flush();f.seek(0);assert hashlib.sha256(f.read()).digest()==hashlib.sha256(raw).digest()
mark('write_flush_reread_hash_existing_paf_seconds',start)
res={'status':'pass','reported':len(records),'duplicate_anchors':0,'independently_scored':scored,'paf_bytes':len(raw),'timings':timings,'caveat':'independent Python audit of emitted sites and cached PAF file copy, not production timings or full native replay'};Path(a.out).write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
