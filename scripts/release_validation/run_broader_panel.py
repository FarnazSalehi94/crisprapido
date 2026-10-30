#!/usr/bin/env python3
"""Index and evaluate a small, already-generated panel using frozen parameters."""
import argparse
from pathlib import Path
import subprocess
import sys
p=argparse.ArgumentParser(description=__doc__)
for name in ('panel','binary','columba','builder','loader','library-path','samtools','matrices','time-bin'):
 p.add_argument('--'+name,required=True)
a=p.parse_args(); d=Path(a.panel).resolve(); ref=str(d/'ref.fa'); guides=str(d/'guides.fa');idx=str(d/'index')
def run(label,cmd):
 with (d/(label+'.stdout')).open('x') as out, (d/(label+'.stderr')).open('x') as err:
  subprocess.run([a.time_bin,'-v','-o',str(d/(label+'.time')),'/usr/bin/timeout','180']+cmd,stdout=out,stderr=err,check=True)
run('build',[a.loader,'--library-path',a.library_path,a.builder,'-f',ref,'-r',idx])
for e,label in ((0,'exact'),(5,'edit')):
 run(label,[a.columba,'-r',idx,'-f',guides,'-a','all','-m','edit','-e',str(e),'-t','1','-o',str(d/(label+'.sam'))])
run('search',[a.binary,'-r',ref,'-G',guides,'-p','GG','-m','4','-b','1','-z','1','-f','0.75','-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-bin',a.columba,'--columba-index',idx])
(d/'out.paf').write_bytes((d/'search.stdout').read_bytes())
subprocess.run([a.samtools,'faidx',ref],check=True)
script=Path(__file__).resolve().parent
with (d/'metrics.jsonl').open('x') as out:subprocess.run([sys.executable,str(script/'compare_indexed_coverage.py'),str(d)],stdout=out,check=True)
subprocess.run([sys.executable,str(script/'resolve_bulge_discrepancies.py'),str(d),'--output',str(d/'alignment_review.json')],check=True)
subprocess.run([sys.executable,str(script/'audit_emitted.py'),'--paf',str(d/'out.paf'),'--fasta',ref,'--guides',guides,'--matrices',a.matrices,'--samtools',a.samtools,'--output',str(d/'audit')],check=True,stdout=subprocess.DEVNULL)
