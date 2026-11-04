#!/usr/bin/env python3
"""Compare integrated PAF representatives to the independently checked native fixtures."""
import csv,subprocess,sys
from pathlib import Path
cli,poc,scratch=sys.argv[1:];p=Path(scratch);root=Path(__file__).resolve().parents[3]
cases=[]
for name in ('indexed_synthetic','indexed_real_chm13','indexed_12_windows'):
 d=root/'docs/coverage_oracles_2026-10-29_v1'/name
 cases.append((d/'ref.fa',d/'guides.fa',p/name,p/(name+'.tsv'),'4','1'))
cases.append((p/'alt/ref.fa',p/'alt/guides.fa',p/'alt/index',p/'alt/out.tsv','0','2'))
for ref,guides,index,tsv,m,z in cases:
 lines=tsv.read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'))
 expected={(r['guide'],r['contig'],r['strand'],r['start'],r['end'],r['cigar']) for r in csv.DictReader(lines[i:],delimiter='\t')}
 cmd=[cli,'--search-backend','pam-anchored','--pam-search-bin',poc,'--columba-index',str(index),'-r',str(ref),'-G',str(guides),'-t','1','-m',m,'-b','1','-z',z]
 out=subprocess.run(cmd,text=True,capture_output=True,check=True,cwd=root).stdout;got=[]
 for line in out.splitlines():
  f=line.split('\t');got.append((f[0],f[5],f[4],f[7],f[8],next(x[5:] for x in f[12:] if x.startswith('cg:Z:'))))
 assert len(got)==len(set(got)) and set(got)==expected,(str(index),len(expected-set(got)),len(set(got)-expected))
 print(str(index),len(got),'representatives matched')
