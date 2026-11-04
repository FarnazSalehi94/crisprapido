#!/usr/bin/env python3
"""Run both complete backends on the same bounded panel and audit each PAF."""
import argparse,json,re,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser()
for x in ('crispr','poc','columba','index','reference','guides','oracle','out'): p.add_argument('--'+x,required=True)
p.add_argument('--repetitions',type=int,default=2);p.add_argument('--time-bin',default='/usr/bin/time')
a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
common=['-r',a.reference,'-G',a.guides,'-m','5','-b','1','-z','2','-f','.75','-p','GG','-t','1','--columba-index',a.index]
root=Path(__file__).resolve().parents[3]
results=[]
for i in range(a.repetitions):
 for backend in ('columba-wfa2','pam-anchored'):
  stem=out/(backend+str(i));cmd=[a.crispr,'--search-backend',backend]+common+(['--columba-bin',a.columba] if backend=='columba-wfa2' else ['--pam-search-bin',a.poc])
  with stem.with_suffix('.paf').open('w') as f,stem.with_suffix('.log').open('w') as e:
   subprocess.run([a.time_bin,'-f','%e %M','-o',str(stem.with_suffix('.time'))]+cmd,stdout=f,stderr=e,check=True,cwd=root)
  audit=subprocess.run(['python3',str(Path(__file__).with_name('audit.py')),a.reference,a.guides,a.oracle,str(stem.with_suffix('.paf'))],text=True,capture_output=True,check=True,cwd=root)
  wall,rss=stem.with_suffix('.time').read_text().split(); log=stem.with_suffix('.log').read_text(); search=re.search(r'states=(\d+) search_wall_seconds=([\d.e+-]+)',log)
  results.append({'backend':backend,'run':i,'wall_seconds':float(wall),'max_process_rss_kib':int(rss),'search_states':int(search[1]) if search else None,'search_seconds':float(search[2]) if search else None,'audit':json.loads(audit.stdout)})
(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
