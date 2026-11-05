#!/usr/bin/env python3
"""Repeat current experimental 2 Mb workflow and compare PAF to saved run."""
import argparse,json,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser()
for x in ('crispr','poc','reference','guides','index','expected','out'):p.add_argument('--'+x,required=True)
p.add_argument('--repetitions',type=int,default=3);a=p.parse_args();root=Path(__file__).resolve().parents[3];expected=Path(a.expected).read_bytes();runs=[]
for i in range(a.repetitions):
 cmd=[a.crispr,'--search-backend','pam-anchored','--pam-search-bin',a.poc,'--columba-index',a.index,'-r',a.reference,'-G',a.guides,'-t','1','-m','5','-b','1','-z','2','-f','.75','-p','GG'];start=time.monotonic();r=subprocess.run(cmd,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90,check=True);elapsed=time.monotonic()-start;assert r.stdout==expected;runs.append({'run':i,'wall_seconds':elapsed,'identical_paf':True})
Path(a.out).write_text(json.dumps(runs,indent=2)+'\n');print(json.dumps(runs))
