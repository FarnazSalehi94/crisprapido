import random
from pathlib import Path
import sys
sys.path.insert(0,'scripts/release_validation')
from coverage_oracle import rc
p=Path(sys.argv[1])
p.mkdir(parents=True, exist_ok=True)
guides={'g1':'ACGTCAGTACGATCGTACGA','g2':'GAGTCCGAGCAGAAGAAGAA'}
r=random.Random(20261030)
def mutate(s,positions):
 s=list(s)
 for i in positions:s[i]={'A':'C','C':'G','G':'T','T':'A'}[s[i]]
 return ''.join(s)
seqs={}
for name,g in guides.items():
 for kind in ['exact','mm1','mm4','I','D','mmI','mmD']:
  t=g
  if kind=='mm1': t=mutate(g,[3])
  if kind=='mm4': t=mutate(g,[2,6,11,16])
  if kind=='I': t=g[:9]+g[10:]
  if kind=='D': t=g[:9]+'T'+g[9:]
  if kind=='mmI': t=mutate(g,[2,6,13,17]); t=t[:9]+t[10:]
  if kind=='mmD': t=mutate(g,[2,6,13,17]); t=t[:9]+'T'+t[9:]
  for strand in ['plus','minus']:
   flank=lambda: ''.join(r.choice('ACGT') for _ in range(60))
   seqs[name+'_'+kind+'_'+strand]=flank()+(t+'AGG' if strand=='plus' else rc(t+'AGG'))+flank()
(p/'ref.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in seqs.items()))
(p/'guides.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in guides.items()))
