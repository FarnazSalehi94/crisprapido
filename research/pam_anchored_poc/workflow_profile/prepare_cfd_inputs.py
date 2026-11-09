#!/usr/bin/env python3
"""Extract emitted targets from FASTA for an isolated production-CFD benchmark."""
import argparse
p=argparse.ArgumentParser();p.add_argument('--paf',required=True);p.add_argument('--reference',required=True);p.add_argument('--fai',required=True);p.add_argument('--out',required=True);a=p.parse_args();rows={}
for line in open(a.fai):
 f=line.rstrip().split('\t');rows[f[0]]=tuple(map(int,f[1:5]))
ref=open(a.reference,'rb');rc=lambda x:x.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def fetch(name,lo,hi):
 length,offset,bases,width=rows[name];assert 0<=lo<=hi<=length;out=[]
 while lo<hi:
  n=min(hi-lo,bases-lo%bases);ref.seek(offset+(lo//bases)*width+lo%bases);out.append(ref.read(n).decode('ascii').upper());lo+=n
 return ''.join(out)
with open(a.out,'w') as out:
 for line in open(a.paf):
  f=line.rstrip().split('\t');name=f[5];lo,hi=int(f[7]),int(f[8]);strand=f[4];seq=fetch(name,lo,hi);seq=seq if strand=='+' else rc(seq);pam=fetch(name,hi,hi+3) if strand=='+' else rc(fetch(name,lo-3,lo));assert pam[0] in 'ACGT' and pam[1:]=='GG';tags={v.split(':',2)[0]:v.split(':',2)[2] for v in f[12:]};out.write('\t'.join([seq,tags['cg'],pam[1:],tags.get('cf','NA')])+'\n')
