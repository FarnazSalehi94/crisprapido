#!/usr/bin/env python3
"""Generate deterministic synthetic and unchanged real-window coverage panels."""
import argparse
import json
from pathlib import Path
import random
import sys
from coverage_oracle import rc

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-windows', required=True)
p.add_argument('--output', required=True)
a=p.parse_args(); out=Path(a.output);out.mkdir(exist_ok=False)

def write_panel(name, refs, guides, provenance):
 d=out/name;d.mkdir()
 (d/'ref.fa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in refs.items()))
 (d/'guides.fa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in guides.items()))
 (d/'panel.json').write_text(json.dumps(provenance,indent=2)+'\n')

def mutate(seq, positions):
 s=list(seq)
 for i in positions:s[i]={'A':'C','C':'G','G':'T','T':'A'}[s[i]]
 return ''.join(s)

guides={'frozen1':'TCTAGCGTTTAACTTCTTCG','frozen2':'TGCCACCAATACCCTTTGGA',
        'repeat':'ACACACACACACACACACAC'}
rng=random.Random(20261031); refs={}; cases=[]
for name,g in guides.items():
 examples=[('exact',g)]
 examples += [('mm_'+str(i),mutate(g,[i])) for i in (0,4,9,14,19)]
 examples += [('four_a',mutate(g,[0,5,10,15])),('four_b',mutate(g,[4,9,14,19]))]
 for i in (0,5,10,15,19):
  examples.append(('I_'+str(i),g[:i]+g[i+1:]))
  examples.append(('D_'+str(i),g[:i]+'T'+g[i:]))
 for i in (5,15):
  m=mutate(g,[1,6,11,18]);examples.append(('mmI_'+str(i),m[:i]+m[i+1:]));examples.append(('mmD_'+str(i),m[:i]+'T'+m[i:]))
 for kind,target in examples:
  for strand in ('plus','minus'):
   flank=lambda: ''.join(rng.choice('ACGT') for _ in range(45))
   site=target+'AGGG'  # overlapping NGG PAMs, retained as separate anchors
   contig=name+'_'+kind+'_'+strand
   refs[contig]=flank()+(site if strand=='plus' else rc(site))+flank()
   cases.append(dict(contig=contig,guide=name,designed=kind,strand=strand,target=target,pam_context='AGGG'))
# Two different exact guides on one overlapping PAM array.
seq='GCTAACGTAGCTAGCATCGT'+'AGGG'
for strand in ('plus','minus'):
 refs['overlap_exact_'+strand]='T'*45+(seq if strand=='plus' else rc(seq))+'C'*45
# A one-base shifted guide preserves exactness for adjacent NGG PAMs.
guides['overlap_a']=seq[:20]; guides['overlap_b']=seq[1:21]
write_panel('synthetic',refs,guides,{'seed':20261031,'cases':cases,'exact_overlap_guides':['overlap_a','overlap_b']})

source=Path(a.source_windows)
source_refs={};name=None
for line in (source/'ref.fa').read_text().splitlines():
 if line.startswith('>'):name=line[1:].split()[0];source_refs[name]=''
 else:source_refs[name]+=line.strip().upper()
selected=['CHM13_window_1','HG04184_mat_window_1','GRCh38_window_1','HG00128_hap2_window_1']
real={k:source_refs[k] for k in selected}
real_guides={}; provenance=[]
for idx,(contig,s) in enumerate(real.items()):
 for i in range(100,len(s)-23):
  if set(s[i:i+23])<=set('ACGT') and s[i+21:i+23]=='GG':break
 else:raise ValueError('no NGG target in '+contig)
 exact=s[i:i+20]
 # 20-nt synthetic query designs compared to unchanged genomic windows.
 variants={'exact':exact,'mm4':mutate(exact,[1,6,11,17]),
           'I':exact[:10]+'A'+exact[10:19],
           'D':s[i-1:i+9]+s[i+10:i+20]}
 for label,guide in variants.items():
  key='r{}_{}'.format(idx+1,label);real_guides[key]=guide
  provenance.append(dict(guide=key,source_contig=contig,source_start=i,design=label,sequence=guide))
write_panel('real',real,real_guides,{'source':'retained 12-assembly windows','designs':provenance,'note':'reference sequences unchanged; guides include deliberate designs'})
