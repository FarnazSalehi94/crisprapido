#!/usr/bin/env python3
"""Compare disposable PAM precheck probe against unchanged CLI across modes."""
import argparse,json,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('panel','output','baseline','probe','columba','matrices'):p.add_argument('--'+n,required=True)
a=p.parse_args();panel=Path(a.panel);out=Path(a.output);out.mkdir(exist_ok=False)
modes=[('m5',['-p','GG','-m','5','-b','1','-z','2','-f','.75']),('m6',['-p','GG','-m','6','-b','1','-z','2','-f','.70']),('m7',['-p','GG','-m','7','-b','1','-z','2','-f','.65']),('strict_ga',['-p','GA','-m','5','-b','1','-z','2','-f','.75']),('no_bulge',['-p','GG','-m','5','-b','0','-z','2','-f','.75']),('size_one',['-p','GG','-m','5','-b','1','-z','1','-f','.75']),('permissive',['-p','GG','-m','5','-b','1','-z','2','-f','.75','--allow-noncanonical-pam']),('no_filter',['-p','GG','-m','5','-b','1','-z','2','-f','.75','--no-filter'])]
results=[]
for name,opts in modes:
 observations=[]
 for label,exe in [('baseline',a.baseline),('probe',a.probe)]:
  cmd=[exe,'-r',str(panel/'ref.fa'),'-G',str(panel/'guides.fa'),'-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-bin',a.columba,'--columba-index',str(panel/'index')]+opts
  t=time.monotonic()
  with (out/f'{name}_{label}.paf').open('w') as o,(out/f'{name}_{label}.stderr').open('w') as e:subprocess.run(cmd,check=True,stdout=o,stderr=e,timeout=120)
  rows=sorted((out/f'{name}_{label}.paf').read_text().splitlines());observations.append((rows,time.monotonic()-t))
 if observations[0][0]!=observations[1][0]:raise ValueError(name+' outputs differ')
 if name in ('m5','m6','m7') and observations[0][0]!=sorted((panel/(name+'.paf')).read_text().splitlines()):raise ValueError(name+' prior oracle-backed records differ')
 row={'mode':name,'rows':len(observations[0][0]),'baseline_seconds':observations[0][1],'probe_seconds':observations[1][1],'records_identical':True};results.append(row);print(json.dumps(row),flush=True)
(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
