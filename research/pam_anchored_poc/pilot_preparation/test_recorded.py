#!/usr/bin/env python3
"""Recheck committed bounded evidence without requiring the large index."""
import json,shutil,subprocess,sys,tempfile
from pathlib import Path
base=Path(__file__).resolve().parent;data=base/'results'
with tempfile.TemporaryDirectory() as d:
 d=Path(d)
 for name in ('stdout.log','stderr.log','measurements.json'):
  source=data/('probe-'+name if name!='measurements.json' else name);shutil.copyfile(source,d/name)
 out=d/'summary.json';subprocess.run([sys.executable,str(base/'check_results.py'),'--cases',str(data/'cases.json'),'--independent',str(data/'independent.json'),'--run',str(d),'--cap','20000','--out',str(out)],check=True,stdout=subprocess.DEVNULL)
 result=json.loads(out.read_text());assert result==json.loads((data/'summary.json').read_text());assert result['selected']==17 and result['completed']==17 and result['skipped']==0
 assert result['completed_by_label_strand']['natural_target_n -']==1
 assert result['completed_by_label_strand']['natural_pam_required_n -']==1
print('17 recorded exact cases and independent source-window checks passed')
