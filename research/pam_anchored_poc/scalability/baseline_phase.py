#!/usr/bin/env python3
"""Separately retain candidate SAM logs to extract Columba mapping-worker elapsed."""
import argparse,hashlib,json,os,re,signal,subprocess,tempfile,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--region',required=True);p.add_argument('--crispr',required=True);p.add_argument('--columba',required=True);p.add_argument('--expected-paf',required=True);p.add_argument('--out',required=True);p.add_argument('--timeout',type=int,default=90);a=p.parse_args();d=Path(a.region)
root=Path(__file__).resolve().parents[3]
with tempfile.TemporaryDirectory(prefix='pam-baseline-phases-') as tmp:
 env=dict(os.environ,TMPDIR=tmp);cmd=[a.crispr,'--search-backend','columba-wfa2','--columba-bin',a.columba,'--columba-index',str(d/'index'),'-r',str(d/'ref.fa'),'-G',str(d/'guides.fa'),'-m','5','-b','1','-z','2','-f','.75','-p','GG','-t','1','--keep-columba-sam']
 start=time.monotonic();proc=subprocess.Popen(['/bin/bash','-c','ulimit -v 2097152; exec "$@"','bounded']+cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=root,env=env,start_new_session=True)
 try:stdout,stderr=proc.communicate(timeout=a.timeout)
 except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.communicate();raise RuntimeError('baseline phase run exceeded limit')
 elapsed=time.monotonic()-start
 if proc.returncode:raise RuntimeError(stderr.decode(errors='replace'))
 expected=Path(a.expected_paf).read_bytes();assert stdout==expected,'instrumented SAM retention changed PAF'
 logs=list(Path(tmp).glob('*/columba.log'));sams=list(Path(tmp).glob('*/output.sam'));times=[]
 for log in logs:times.extend(float(x) for x in re.findall(r'mapping threads finished mapping in ([0-9.]+)s',log.read_text()))
 result={'wall_seconds':elapsed,'columba_processes':len(logs),'mapping_worker_seconds':times,'mapping_worker_seconds_sum':sum(times) if len(times)==len(logs) and logs else None,'sam_bytes':sum(x.stat().st_size for x in sams),'paf_sha256':hashlib.sha256(stdout).hexdigest(),'matches_uninstrumented_paf':True,'limitations':'worker time includes mapper reporting and is not equivalent to the prototype timer; SAM write time is not isolated'}
 Path(a.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
