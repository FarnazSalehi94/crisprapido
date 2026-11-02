#!/usr/bin/env python3
"""Submit four serialized low-priority jobs after successful verification."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',required=True);p.add_argument('--verification-job',required=True,type=int)
a=p.parse_args();root=Path(a.root).resolve()
if (root/'submitted.json').exists(): raise SystemExit('already submitted')
rows={x['assembly_id']:x for x in csv.DictReader((root/'inputs/verified_assemblies.tsv').open(),delimiter='\t')}
previous=a.verification_job;submitted=[]
try:
 for assembly in ('CHM13','GRCh38'):
  for guide in ('ELMOD2_TCTAGCGTTTAACTTCTTCG','MED6_TGCCACCAATACCCTTTGGA'):
   r=rows[assembly]
   env=dict(os.environ,PILOT_ROOT=str(root),PILOT_ASSEMBLY=assembly,PILOT_GUIDE=guide,PILOT_REFERENCE=r['fasta'],PILOT_INDEX=r['index_prefix'])
   cmd=['sbatch','--parsable',f'--dependency=afterok:{previous}',f'--output={root}/search.%j.out',f'--error={root}/search.%j.err',str(root/'tools/slurm_whole_assembly_pilot.sbatch')]
   job=int(subprocess.check_output(cmd,env=env,text=True).strip().split(';')[0])
   submitted.append({'job':job,'afterok':previous,'assembly':assembly,'guide':guide});previous=job
   (root/'submitted.partial.json').write_text(json.dumps(submitted,indent=2)+'\n')
except Exception:
 if submitted: subprocess.run(['scancel']+[str(x['job']) for x in submitted],check=False)
 raise
(root/'submitted.json').write_text(json.dumps(submitted,indent=2)+'\n')
print(json.dumps(submitted,indent=2))
