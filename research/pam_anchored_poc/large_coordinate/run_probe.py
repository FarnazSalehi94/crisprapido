#!/usr/bin/env python3
"""Run bounded probe with sampled process-group memory and fail-closed limits."""
import argparse,json,os,signal,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--probe',required=True);p.add_argument('--stage',required=True);p.add_argument('--cases',required=True);p.add_argument('--out',required=True);p.add_argument('--timeout',type=int,default=900);p.add_argument('--aggregate-kib',type=int,default=23068672);p.add_argument('--vm-kib',type=int,default=25165824);a=p.parse_args();stage=Path(a.stage);out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
def group_rss(group):
 total=0;members=0
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   if os.getpgid(int(p.name))!=group:continue
   for line in (p/'status').read_text().splitlines():
    if line.startswith('VmRSS:'):total+=int(line.split()[1]);members+=1;break
  except (ProcessLookupError,PermissionError,FileNotFoundError):pass
 return total,members
cmd=['/bin/bash','-c','ulimit -v "$1"; shift; exec "$@"','limited',str(a.vm_kib),a.probe,str(stage/'index'),str(stage/'reference.fa'),a.cases];start=time.monotonic();peak=0;phase='startup';peaks={};events={};reason=None;seen=0
with (out/'stdout.log').open('w') as fo,(out/'stderr.log').open('w') as fe:
 proc=subprocess.Popen(cmd,stdout=fo,stderr=fe,start_new_session=True)
 with (out/'stderr.log').open('r') as monitor:
  while proc.poll() is None:
   for line in monitor:
    fields=line.rstrip().split('\t')
    if len(fields)==2 and fields[0]=='PHASE':phase=fields[1]
    elif len(fields)==3 and fields[0]=='EVENT':events[fields[1]]=float(fields[2])
   rss,n=group_rss(proc.pid);peak=max(peak,rss);peaks[phase]=max(peaks.get(phase,0),rss);seen+=1
   if rss>a.aggregate_kib or time.monotonic()-start>a.timeout:
    reason='aggregate_rss_limit' if rss>a.aggregate_kib else 'timeout';os.killpg(proc.pid,signal.SIGKILL);break
   time.sleep(.1)
  proc.wait()
  for line in monitor:
   f=line.rstrip().split('\t')
   if len(f)==3 and f[0]=='EVENT':events[f[1]]=float(f[2])
elapsed=time.monotonic()-start;result={'returncode':proc.returncode,'stop_reason':reason,'wall_seconds':elapsed,'sample_interval_seconds':.1,'sample_count':seen,'sampled_process_group_rss_peak_kib':peak,'sampled_phase_peaks_kib':peaks,'events':events,'memory_note':'sum of process RSS; can double-count shared pages and miss peaks between samples','wall_limit_seconds':a.timeout,'aggregate_rss_limit_kib':a.aggregate_kib,'per_process_vm_limit_kib':a.vm_kib}
(out/'measurements.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));
if proc.returncode:raise SystemExit(1)
