#!/usr/bin/env python3
"""Bounded repeated integrated backend comparison; no indexes or references committed."""
import argparse,hashlib,json,os,re,signal,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser()
for n in ('dataset','crispr','poc','columba','time_bin','out'):p.add_argument('--'+n.replace('_','-'),required=True)
p.add_argument('--repetitions',type=int,default=3);p.add_argument('--timeout',type=int,default=90);p.add_argument('--memory-kib',type=int,default=2097152)
a=p.parse_args();root=Path(__file__).resolve().parents[3];dataset=Path(a.dataset);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
regions=json.loads((dataset/'regions.json').read_text());results=[]
for region in regions:
 name=region['name'];d=dataset/name;dest=out/name;dest.mkdir(exist_ok=True)
 for rep in range(a.repetitions):
  for backend in ('columba-wfa2','pam-anchored'):
   stem=dest/(backend+'-'+str(rep));common=['-r',str(d/'ref.fa'),'-G',str(d/'guides.fa'),'-m','5','-b','1','-z','2','-f','0.75','-p','GG','-t','1','--columba-index',str(d/'index')]
   cmd=[a.crispr,'--search-backend',backend]+common+(['--columba-bin',a.columba] if backend=='columba-wfa2' else ['--pam-search-bin',a.poc])
   bounded=['/bin/bash','-c','ulimit -v "$1"; shift; exec "$@"','bounded',str(a.memory_kib),a.time_bin,'-f','%e %M','-o',str(stem)+'.time']+cmd
   with open(str(stem)+'.paf','w') as fo,open(str(stem)+'.log','w') as fe:
    proc=subprocess.Popen(bounded,cwd=root,stdout=fo,stderr=fe,start_new_session=True)
    try:code=proc.wait(timeout=a.timeout)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait();code='timeout'
   log=Path(str(stem)+'.log').read_text(); t=Path(str(stem)+'.time'); fields=t.read_text().strip().splitlines()[-1].split() if t.exists() and t.stat().st_size else []
   search=re.search(r'states=(\d+) search_wall_seconds=([\deE.+-]+)',log)
   result={'region':name,'run':rep,'backend':backend,'status':code,'wall_seconds':float(fields[0]) if len(fields)==2 else None,'max_individual_process_rss_kib':int(fields[1]) if len(fields)==2 else None,'states':int(search[1]) if search else None,'search_seconds':float(search[2]) if search else None,'rows':sum(1 for _ in open(str(stem)+'.paf')) if code==0 else None}
   results.append(result);(out/'runs.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(result),flush=True)
   if code!=0:
    print('STOP: resource limit or backend failure; no further repetitions for this region',flush=True);break
  if code!=0:break
