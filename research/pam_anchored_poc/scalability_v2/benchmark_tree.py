#!/usr/bin/env python3
"""Compare bounded workflows, sampling process-group RSS and enforcing limits."""
import argparse,json,os,re,signal,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser()
for x in ('dataset','crispr','poc','columba','out'):p.add_argument('--'+x,required=True)
p.add_argument('--timeout',type=int,default=90);p.add_argument('--aggregate-kib',type=int,default=1572864);p.add_argument('--repetitions',type=int,default=1);a=p.parse_args();dataset=Path(a.dataset);out=Path(a.out);out.mkdir(parents=True,exist_ok=True);root=Path(__file__).resolve().parents[3];records=[]
def group_rss(group):
 total=0;count=0
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   if os.getpgid(int(p.name))!=group:continue
   for line in (p/'status').read_text().splitlines():
    if line.startswith('VmRSS:'):total+=int(line.split()[1]);count+=1;break
  except (ProcessLookupError,PermissionError,FileNotFoundError):pass
 return total,count
for d in (dataset/'CHM13',dataset/'GRCh38'):
 for rep in range(a.repetitions):
  for backend in ('columba-wfa2','pam-anchored'):
   stem=out/(d.name+'-'+backend+'-'+str(rep));cmd=[a.crispr,'--search-backend',backend,'-r',str(d/'ref.fa'),'-G',str(d/'guides.fa'),'-m','5','-b','1','-z','2','-f','.75','-p','GG','-t','1','--columba-index',str(d/'index')]+(['--columba-bin',a.columba] if backend=='columba-wfa2' else ['--pam-search-bin',a.poc]);cmd=['/bin/bash','-c','ulimit -v 2097152; exec "$@"','bounded']+cmd
   start=time.monotonic();peak=0;peak_processes=0;reason=None
   with open(str(stem)+'.paf','w') as fo,open(str(stem)+'.log','w') as fe:
    proc=subprocess.Popen(cmd,stdout=fo,stderr=fe,cwd=root,start_new_session=True)
    while proc.poll() is None:
     rss,n=group_rss(proc.pid)
     if rss>peak:peak=rss;peak_processes=n
     if rss>a.aggregate_kib or time.monotonic()-start>a.timeout:
      reason='aggregate_rss_limit' if rss>a.aggregate_kib else 'timeout';os.killpg(proc.pid,signal.SIGKILL);break
     time.sleep(.1)
    proc.wait()
   elapsed=time.monotonic()-start;log=Path(str(stem)+'.log').read_text();m=re.search(r'states=(\d+) search_wall_seconds=([\deE.+-]+)',log)
   item={'region':d.name,'backend':backend,'run':rep,'status':proc.returncode,'stop_reason':reason,'wall_seconds':elapsed,'sampled_aggregate_rss_kib':peak,'sampled_processes_at_peak':peak_processes,'states':int(m[1]) if m else None,'search_seconds':float(m[2]) if m else None,'rows':sum(1 for _ in open(str(stem)+'.paf')) if proc.returncode==0 else None};records.append(item);(out/'runs.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(item),flush=True)
   if proc.returncode:raise SystemExit('resource limit or backend failure: stopping study')
