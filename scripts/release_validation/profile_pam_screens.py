#!/usr/bin/env python3
"""Count WFA calls using a disposable instrumented build and verify output parity."""
import argparse,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for x in ('benchmark','panel','profile','matrices','output'):p.add_argument('--'+x,required=True)
a=p.parse_args();b=Path(a.benchmark);panel=Path(a.panel);out=Path(a.output);out.mkdir(exist_ok=False);result={}
for mode in ('baseline','overlap_pam'):
 cmd=[a.profile,'-r',str(panel/'ref.fa'),'-G',str(panel/'guides.fa'),'-p','GG','-m','5','-b','1','-z','2','-f','.75','-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-sam',str(b/(mode+'.sam'))]
 with (out/(mode+'.paf')).open('x') as o,(out/(mode+'.stderr')).open('x') as e:subprocess.run(cmd,check=True,stdout=o,stderr=e,timeout=300)
 expected=sorted((b/(mode+'.paf')).read_text().splitlines());actual=sorted((out/(mode+'.paf')).read_text().splitlines())
 if expected!=actual:raise ValueError('instrumented PAF differs in '+mode)
 profile={}
 for line in (out/(mode+'.stderr')).read_text().splitlines():
  f=line.split('\t')
  if len(f)==4 and f[0]=='PROFILE':profile[f[1]]={'calls':int(f[2]),'nanoseconds':int(f[3])}
 if 'wfa_align' not in profile:raise ValueError('missing WFA counter')
 result[mode]=profile
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({x:result[x]['wfa_align'] for x in result}))
