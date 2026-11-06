#!/usr/bin/env python3
"""Test process-tree sampling and fail-closed memory cap without an FM index."""
import json,subprocess,sys,tempfile
from pathlib import Path
here=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='probe-monitor-test-') as td:
 p=Path(td);probe=p/'probe.py';probe.write_text('#!/usr/bin/python3\nimport sys,time\nsys.stderr.write("PHASE\\tindex\\n");sys.stderr.flush();time.sleep(.3)\nsys.stderr.write("PHASE\\tlookup\\n");sys.stderr.flush();time.sleep(.3)\n');probe.chmod(0o755)
 def run(name,limit):
  return subprocess.run([sys.executable,str(here/'run_probe.py'),'--probe',str(probe),'--stage',str(p),'--cases',str(p/'x'),'--out',str(p/name),'--timeout','10','--aggregate-kib',str(limit),'--vm-kib','1048576'],capture_output=True,text=True,timeout=20)
 a=run('ok',1000000);assert a.returncode==0,a.stderr;x=json.loads((p/'ok/measurements.json').read_text());assert x['sampled_process_group_rss_peak_kib']>0 and x['sampled_phase_peaks_kib'].get('index',0)>0 and x['sampled_phase_peaks_kib'].get('lookup',0)>0
 b=run('stop',1);assert b.returncode!=0;y=json.loads((p/'stop/measurements.json').read_text());assert y['stop_reason']=='aggregate_rss_limit'
print('process-group sampling and memory-limit stop passed')
