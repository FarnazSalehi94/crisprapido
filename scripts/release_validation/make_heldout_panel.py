#!/usr/bin/env python3
"""Generate a held-out synthetic panel and four fresh genomic windows."""
import argparse
import json
from pathlib import Path
import random
import subprocess
from coverage_oracle import rc
p=argparse.ArgumentParser(description=__doc__)
for name in ('output','chm13','grch38','samtools'):p.add_argument('--'+name,required=True)
a=p.parse_args();root=Path(a.output);root.mkdir(exist_ok=False)
def write(name,refs,guides,details):
 d=root/name;d.mkdir();(d/'ref.fa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in refs.items()))
 (d/'guides.fa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in guides.items()))
 (d/'panel.json').write_text(json.dumps(details,indent=2)+'\n')
def mut(s,ii):
 x=list(s)
 for i in ii:x[i]={'A':'G','C':'T','G':'C','T':'A'}[x[i]]
 return ''.join(x)
r=random.Random(20261102)
guides={'held_a':'AAATTCACTTACCTAAAACC','held_b':'TCTGCATCAGAGTCATATGC','held_repeat':'GAGAGAGAGAGAGAGAGAGA'};refs={};cases=[]
for name,g in guides.items():
 variants=[('exact',g)] + [('mm'+str(i),mut(g,[i])) for i in (1,3,8,13,18)]
 variants += [('mm4',mut(g,[1,5,12,18]))]
 for i in (1,4,9,14,18):
  variants += [('I'+str(i),g[:i]+g[i+1:]),('D'+str(i),g[:i]+'C'+g[i:])]
 m=mut(g,[2,6,12,17])
 for i in (4,14):variants += [('mixI'+str(i),m[:i]+m[i+1:]),('mixD'+str(i),m[:i]+'C'+m[i:])]
 for typ,seq in variants:
  for strand in ('plus','minus'):
   flank=lambda:''.join(r.choice('ACGT') for _ in range(70))
   core=seq+'TGGG';key=name+'_'+typ+'_'+strand
   refs[key]=flank()+(core if strand=='plus' else rc(core))+flank()
   cases.append({'contig':key,'guide':name,'design':typ,'strand':strand})
write('synthetic',refs,guides,{'seed':20261102,'cases':cases})
refs={};metadata=[]
for label,path in [('CHM13',a.chm13),('GRCh38',a.grch38)]:
 for j,line in enumerate(Path(path+'.fai').read_text().splitlines()[:2]):
  chrom,length=line.split('\t')[:2];first=500001+400000*j;last=first+999
  if int(length)<last:raise ValueError('short contig')
  region=f'{chrom}:{first}-{last}'
  data=subprocess.run([a.samtools,'faidx',path,region],stdout=subprocess.PIPE,universal_newlines=True,check=True).stdout
  seq=''.join(data.splitlines()[1:]).upper()
  if len(seq)!=1000:raise ValueError('truncated region')
  key=f'{label}_new_{j+1}';refs[key]=seq;metadata.append({'name':key,'fasta':path,'region':region})
guides={};designs=[]
for n,(chrom,seq) in enumerate(refs.items()):
 for i in range(100,950):
  t=seq[i:i+23]
  if len(t)==23 and set(t)<=set('ACGT') and t[-2:]=='GG':break
 else:raise ValueError('no PAM')
 g=t[:20]; names={'exact':g,'mm4':mut(g,[0,7,12,18]),'I':g[:9]+'T'+g[9:19],'D':seq[i-1:i+8]+seq[i+9:i+20]}
 for kind,s in names.items():
  key=f'held_real{n}_{kind}';guides[key]=s;designs.append({'guide':key,'contig':chrom,'offset':i,'type':kind})
write('real',refs,guides,{'windows':metadata,'designs':designs})
