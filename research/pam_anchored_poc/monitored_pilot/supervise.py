#!/usr/bin/env python3
"""Supervise a process tree; sampled RSS is not a kernel-enforced limit."""
import argparse, ctypes, json, os, signal, subprocess, sys, time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--limit-kib',type=int,required=True);p.add_argument('--interval',type=float,default=.1);p.add_argument('--log',required=True);p.add_argument('--summary',required=True);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
if a.limit_kib<=0 or not 0<a.interval<=1 or not a.command:raise SystemExit('invalid supervisor arguments')
if a.command[0]=='--':a.command=a.command[1:]
# Adopt orphan descendants: all remain visible to the process-tree monitor.
if ctypes.CDLL(None,use_errno=True).prctl(36,1,0,0,0)!=0:raise SystemExit('cannot become child subreaper')
stop_signal=None
def on_signal(sig,frame):
 global stop_signal
 stop_signal=sig
for sig in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(sig,on_signal)
root=os.getpid();proc=None;peak=0;count=0;reason=None;start=time.monotonic();last_records={};last_log=0

def snapshot():
 records={}
 for path in Path('/proc').iterdir():
  if not path.name.isdigit():continue
  try:
   lines=(path/'status').read_text().splitlines();fields={x.split(':',1)[0]:x.split(':',1)[1].strip() for x in lines if ':' in x}
   records[int(path.name)]=(int(fields['PPid']),int(fields.get('VmRSS','0 kB').split()[0]))
  except (FileNotFoundError,ProcessLookupError):continue
  except PermissionError:continue
 if root not in records:raise RuntimeError('cannot observe supervisor in /proc')
 members={root}
 while True:
  new=members|{pid for pid,(ppid,rss) in records.items() if ppid in members}
  if new==members:break
  members=new
 return records,members,sum(records[pid][1] for pid in members)

def kill_tree():
 # SIGSTOP first limits racing forks; SIGKILL also terminates stopped processes.
 for _ in range(4):
  try: records,members,_=snapshot()
  except Exception: members=set(last_records)
  members.discard(root)
  for sig in (signal.SIGSTOP,signal.SIGKILL):
   if proc:
    try:os.killpg(proc.pid,sig)
    except ProcessLookupError:pass
   for pid in members:
    try:os.kill(pid,sig)
    except ProcessLookupError:pass
  time.sleep(.05)
try:
 with open(a.log,'x',buffering=1) as log:
  log.write('timestamp elapsed_seconds members summed_rss_kib sampled_peak_kib\n')
  proc=subprocess.Popen(a.command,start_new_session=True)
  while True:
   records,members,rss=snapshot();last_records={pid:records[pid] for pid in members if pid!=root};peak=max(peak,rss);count+=1;now=time.monotonic()
   if now-last_log>=1 or rss>a.limit_kib:
    log.write('%s %.3f %d %d %d\n'%(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),now-start,len(members),rss,peak));last_log=now
   if rss>a.limit_kib:reason='sampled_rss_limit';break
   if stop_signal:reason='signal_%d'%stop_signal;break
   code=proc.poll()
   if code is not None:
    # A process that leaves a child behind cannot be considered complete.
    records,members,_=snapshot()
    if len(members)>1:reason='orphan_descendant';break
    reason='completed' if code==0 else 'child_exit_%d'%code;break
   time.sleep(a.interval)
except BaseException as exc:
 reason='monitor_error: '+str(exc)
finally:
 if proc and reason!='completed':kill_tree()
 if proc:
  try:proc.wait(timeout=5)
  except subprocess.TimeoutExpired:pass
 summary={'reason':reason,'sampled_peak_rss_kib':peak,'limit_kib':a.limit_kib,'sample_count':count,'sample_interval_seconds':a.interval,'elapsed_seconds':time.monotonic()-start,'child_returncode':None if not proc else proc.returncode,'caveat':'sampled summed process-tree RSS; may miss spikes and double-count shared pages; not a hard memory limit'}
 try:Path(a.summary).write_text(json.dumps(summary,indent=2)+'\n')
 except OSError as e:print('cannot write supervisor summary: '+str(e),file=sys.stderr);sys.exit(1)
 print(json.dumps(summary),flush=True)
if reason!='completed':sys.exit(1)
