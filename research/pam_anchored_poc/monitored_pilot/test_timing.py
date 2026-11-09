#!/usr/bin/env python3
"""Exercise Bash timing and error propagation without an index or guide search."""
import subprocess,tempfile
from pathlib import Path
helper=Path(__file__).with_name('portable_timing.sh').resolve()
with tempfile.TemporaryDirectory() as td:
 td=Path(td)
 def run(name, command):
  result=td/(name+'.time')
  script='set -Eeuo pipefail; source "$1"; SECONDS=0; run_timed "$2" bash -c "$3"; printf after'
  p=subprocess.run(['bash','-c',script,'timing-test',str(helper),str(result),command],capture_output=True,text=True,timeout=10)
  fields=dict(x.split('=',1) for x in result.read_text().splitlines())
  return p,fields
 p,f=run('success','sleep 1');assert p.returncode==0 and p.stdout=='after' and int(f['elapsed_seconds'])>=1 and f['exit_code']=='0'
 p,f=run('failure','exit 37');assert p.returncode==37 and p.stdout=='' and f['exit_code']=='37' and int(f['elapsed_seconds'])>=0
 p,f=run('missing','exec /definitely/not/a/binary');assert p.returncode==127 and f['exit_code']=='127'
print('Bash timing records success, failure and missing-command status; failures stop the workload')
