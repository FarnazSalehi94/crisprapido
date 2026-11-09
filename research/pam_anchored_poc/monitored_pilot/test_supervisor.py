#!/usr/bin/env python3
import json,subprocess,sys,tempfile
from pathlib import Path
s=Path(__file__).with_name('supervise.py')
with tempfile.TemporaryDirectory() as d:
 d=Path(d)
 def run(name,limit,cmd):
  log=d/(name+'.log');summary=d/(name+'.json');r=subprocess.run([sys.executable,str(s),'--limit-kib',str(limit),'--interval','.05','--log',str(log),'--summary',str(summary),'--']+cmd,capture_output=True,text=True,timeout=15);return r,json.loads(summary.read_text())
 r,x=run('ok',1024*1024,[sys.executable,'-c','print("ok")']);assert r.returncode==0 and x['reason']=='completed',(r.stdout,r.stderr,x)
 r,x=run('fail',1024*1024,[sys.executable,'-c','raise SystemExit(4)']);assert r.returncode!=0 and x['reason']=='child_exit_4',(r.stdout,r.stderr,x)
 r,x=run('limit',30000,[sys.executable,'-c','import subprocess,sys,time;p=subprocess.Popen([sys.executable,"-c","import time; x=bytearray(80*1024*1024);time.sleep(10)"]);time.sleep(10)']);assert r.returncode!=0 and x['reason']=='sampled_rss_limit' and x['sampled_peak_rss_kib']>30000,(r.stdout,r.stderr,x)
print('supervisor normal exit, child failure, and process-tree limit termination passed')
