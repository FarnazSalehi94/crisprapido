#!/usr/bin/env python3
"""Build bounded 2 Mb regional indexes across two assemblies and two strata."""
import argparse
import json
from pathlib import Path
import resource
import subprocess
p=argparse.ArgumentParser(description=__doc__)
for n in ('output','chm13','grch38','guides','samtools','builder','loader','library-path'):p.add_argument('--'+n,required=True)
a=p.parse_args();root=Path(a.output);root.mkdir(exist_ok=False)
# Same explicit, easily audited positions on each assembly; no genome-wide search.
windows=[('chr1','100000001-102000000','nonrepeat_candidate'),('chr9','34000001-36000000','repeat_candidate')]
rows=[]
for assembly,fasta in [('CHM13',a.chm13),('GRCh38',a.grch38)]:
 contigs={line.split('\t')[0].split('#')[-1]:line.split('\t')[0] for line in Path(fasta+'.fai').read_text().splitlines()}
 for chrom,span,stratum in windows:
  name=assembly+'_'+stratum;d=root/name;d.mkdir();region=contigs[chrom]+':'+span
  with (d/'ref.fa').open('x') as out:
   subprocess.run([a.samtools,'faidx',fasta,region],stdout=out,check=True,timeout=60)
  seq=''.join((d/'ref.fa').read_text().splitlines()[1:]);seq_upper=seq.upper()
  if len(seq)!=2000000:raise ValueError('incomplete window '+region)
  acgt=sum(c in 'ACGT' for c in seq_upper);n=sum(c not in 'ACGT' for c in seq_upper)
  di=sum(seq_upper[i]==seq_upper[i-2] and seq_upper[i] in 'ACGT' for i in range(2,len(seq_upper)))
  with (d/'build.log').open('x') as log:
   subprocess.run([a.loader,'--library-path',a.library_path,a.builder,'-f',str(d/'ref.fa'),'-r',str(d/'index')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
  rows.append({'name':name,'assembly':assembly,'region':region,'stratum_label':stratum,'bases':len(seq),'acgt':acgt,'ambiguous':n,'period2_matches':di,'period2_denominator':len(seq)-2})
(root/'windows.json').write_text(json.dumps(rows,indent=2)+'\n')
lines=Path(a.guides).read_text().splitlines();(root/'guides.fa').write_text('\n'.join(lines[:8])+'\n')
