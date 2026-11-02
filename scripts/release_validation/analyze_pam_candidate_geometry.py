#!/usr/bin/env python3
"""Exact external-sort SAM key analysis and strict overlapping-PAM feasibility.

Reads an existing SAM; launches no search. Sorting scratch files are removed on
success. This analysis is specific to 20bp guides, strict NGG and b=1/z=2.
"""
import argparse, collections, json, os, re, subprocess, time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for x in ('sam','exact-sam','fasta','samtools','output'):p.add_argument('--'+x,required=True)
p.add_argument('--sort-memory-mib',type=int,default=512)
a=p.parse_args();t0=time.monotonic();out=Path(a.output);out.mkdir(exist_ok=False)
if a.sort_memory_mib<1:raise ValueError('invalid sort memory')
exact=set()
with open(a.exact_sam) as f:
 for line in f:
  if line.startswith('@'):continue
  x=line.split('\t');flag=int(x[1])
  if not flag&4:exact.add((x[2],x[0],int(x[3])-1,bool(flag&16),x[5]))
keys=out/'keys.unsorted';sortedkeys=out/'keys.sorted';nm=collections.Counter();unmapped=raw_count=0
with open(a.sam) as src,keys.open('w') as dst:
 for line in src:
  if line.startswith('@'):continue
  x=line.rstrip('\n').split('\t');flag=int(x[1])
  if flag&4:unmapped+=1;continue
  nmtag=next((v for v in x[11:] if v.startswith('NM:i:')),None)
  if nmtag is None:raise ValueError('missing NM tag')
  nm[nmtag[5:]]+=1;raw_count+=1
  dst.write(f'{x[2]}\t{x[0]}\t{int(x[3])-1}\t{int(bool(flag&16))}\t{x[5]}\n')
  if raw_count%1000000==0:print('extracted',raw_count,flush=True)
print('sorting',raw_count,flush=True)
env=dict(os.environ,LC_ALL='C')
subprocess.run(['sort','--parallel=1','-S',str(a.sort_memory_mib)+'M','-T',str(out),'-o',str(sortedkeys),str(keys)],check=True,env=env)
totals=collections.Counter();unique=collections.Counter();found_exact=set();seq=b'';last_contig=None;DNA=b'ACGT'
def getseq(contig):
 result=subprocess.run([a.samtools,'faidx',a.fasta,contig],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 return b''.join(result.stdout.splitlines()[1:]).upper()
def valid_in_range(seq,lo,hi,rev):
 if lo>hi:return False
 pair=b'CC' if rev else b'GG';cursor=lo if rev else lo+1;limit=hi+2 if rev else hi+3
 while True:
  i=seq.find(pair,cursor,limit)
  if i<0:return False
  p=i if rev else i-1
  if lo<=p<=hi and seq[p+2 if rev else p] in DNA:return True
  cursor=i+1
def process(line,n):
 global last_contig,seq
 x=line.rstrip('\n').split('\t');contig,name,pos,rev,cigar=x[0],x[1],int(x[2]),bool(int(x[3])),x[4]
 if contig!=last_contig:seq=getseq(contig);last_contig=contig;print('reference',contig,flush=True)
 refspan=sum(int(v) for v,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar) if op in 'MDN=X');end=pos+refspan
 if refspan==0 or end>len(seq):raise ValueError('invalid SAM reference span')
 if rev:
  direct=pos>=3 and seq[pos-3:pos-1]==b'CC' and seq[pos-1] in DNA;lo=max(0,pos-24);hi=min(end-4,len(seq)-21)
 else:
  direct=end+3<=len(seq) and seq[end+1:end+3]==b'GG' and seq[end] in DNA;lo=max(18,pos+1);hi=min(end+21,len(seq)-3)
 possible=valid_in_range(seq,lo,hi,rev)
 if direct and 18<=refspan<=22 and not possible:raise ValueError('internal geometry contradiction')
 fullkey=(contig,name,pos,rev,cigar)
 if fullkey in exact:found_exact.add(fullkey)
 for counts,weight in ((totals,n),(unique,1)):
  counts['mapped']+=weight;counts['direct_pam_valid']+=weight*direct;counts['direct_pam_invalid']+=weight*(not direct);counts['no_possible_overlapping_pam']+=weight*(not possible)
previous=None;count=0
with sortedkeys.open() as f:
 for line in f:
  if line!=previous:
   if previous is not None:process(previous,count)
   previous=line;count=1
  else:count+=1
 if previous is not None:process(previous,count)
if totals['mapped']!=raw_count:raise ValueError('raw count mismatch after sort')
summary={'sam':str(Path(a.sam).resolve()),'mapped_raw':totals['mapped'],'mapped_unique':unique['mapped'],'duplicate_keys':totals['mapped']-unique['mapped'],'exact_pass_keys':len(exact),'exact_pass_keys_already_in_e7':len(found_exact),'e7_nm_distribution':dict(sorted(nm.items())),'raw':dict(totals),'unique':dict(unique),'unmapped':unmapped,'elapsed_seconds':time.monotonic()-t0,'sort_memory_mib':a.sort_memory_mib,'definition':'direct PAM is adjacent to the raw candidate span; no-overlapping-PAM tests all 18..22bp spans overlapping it on the same strand, with complete strict NGG; no other candidate rejection is assumed'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');keys.unlink();sortedkeys.unlink();print(json.dumps(summary,sort_keys=True))
