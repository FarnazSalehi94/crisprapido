#!/usr/bin/env python3
"""Compare old and development executables on the same small indexed fixtures."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess
import time
p = argparse.ArgumentParser(description=__doc__)
for name in ('panels', 'old', 'new', 'columba', 'output'):
    p.add_argument('--'+name, required=True)
p.add_argument('--repeats', type=int, default=10)
a = p.parse_args()
if a.repeats < 1: p.error('repeats must be positive')
root = Path(a.panels).resolve()
results=[]
for panel in ('indexed_synthetic', 'indexed_real_chm13', 'indexed_12_windows'):
    d=root/panel; samples={'old':[], 'new':[]}; counts={}
    for i in range(a.repeats):
        # Alternate order to reduce systematic warm-cache bias.
        for label in (('old','new') if i%2 == 0 else ('new','old')):
            binary=getattr(a,label)
            cmd=[binary,'-r',str(d/'ref.fa'),'-G',str(d/'guides.fa'),'-p','GG','-m','4','-b','1','-z','1','-f','0.75','-t','1','--columba-bin',a.columba,'--columba-index',str(d/'index')]
            start=time.monotonic()
            result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
            samples[label].append(time.monotonic()-start)
            count=len(result.stdout.splitlines())
            if label in counts and counts[label]!=count: raise ValueError('nondeterministic row count')
            counts[label]=count
    medians={k:statistics.median(v) for k,v in samples.items()}
    results.append({'panel':panel,'repeats':a.repeats,'seconds':samples,'median_seconds':medians,
                    'new_over_old':medians['new']/medians['old'],'row_counts':counts})
Path(a.output).write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results))
