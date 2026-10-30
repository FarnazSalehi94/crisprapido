#!/usr/bin/env python3
"""Alternating bounded indexed-pilot runs; retain PAF, stderr and time per run."""
import argparse
import json
from pathlib import Path
import subprocess
import time
p=argparse.ArgumentParser(description=__doc__)
for n in ('panel','output','base','fast','columba','time-bin'):p.add_argument('--'+n,required=True)
p.add_argument('--repeats',type=int,default=2);a=p.parse_args()
if a.repeats<1:p.error('positive repeats required')
d=Path(a.panel).resolve();out=Path(a.output);out.mkdir(exist_ok=False);samples=[];expected=None
for i in range(a.repeats):
 for label in (('base','fast') if i%2==0 else ('fast','base')):
  prefix=out/('{}_{}'.format(label,i));paf=prefix.with_suffix('.paf');err=prefix.with_suffix('.stderr');tim=prefix.with_suffix('.time')
  cmd=[a.time_bin,'-v','-o',str(tim),'/usr/bin/timeout','180',getattr(a,label),'-r',str(d/'ref.fa'),'-G',str(d/'guides.fa'),'-p','GG','-m','4','-b','1','-z','1','-f','0.75','-t','1','--columba-bin',a.columba,'--columba-index',str(d/'index')]
  start=time.monotonic()
  with paf.open('xb') as stdout,err.open('xb') as stderr:subprocess.run(cmd,check=True,stdout=stdout,stderr=stderr)
  elapsed=time.monotonic()-start;data=paf.read_bytes()
  if expected is None:expected=data
  elif data!=expected:raise ValueError('output changed between benchmark runs')
  samples.append({'iteration':i,'binary':label,'elapsed_seconds':elapsed,'paf_rows':len(data.splitlines()),'time_file':str(tim)})
(out/'samples.json').write_text(json.dumps(samples,indent=2)+'\n')
