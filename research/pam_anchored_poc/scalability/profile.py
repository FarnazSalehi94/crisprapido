#!/usr/bin/env python3
"""Measure isolated copy/check costs and native startup overhead on bounded input."""
import argparse,json,shutil,subprocess,sys,tempfile,time,re
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'compatibility'))
from index_manifest import SUFFIXES,verify
p=argparse.ArgumentParser();p.add_argument('--dataset',required=True);p.add_argument('--poc',required=True);p.add_argument('--out',required=True);p.add_argument('--repetitions',type=int,default=3);p.add_argument('--region',action='append');a=p.parse_args();dataset=Path(a.dataset);results=[]
for region in json.loads((dataset/'regions.json').read_text()):
 if a.region and region['name'] not in a.region:continue
 d=dataset/region['name']
 for rep in range(a.repetitions):
  with tempfile.TemporaryDirectory(prefix='pam-scaling-profile-') as td:
   t=Path(td);start=time.monotonic();shutil.copyfile(d/'ref.fa',t/'ref.fa')
   for suffix in (*SUFFIXES,'crispr-binding.json'):shutil.copyfile(str(d/'index')+'.'+suffix,str(t/'index')+'.'+suffix)
   copied=time.monotonic();verify(t/'ref.fa',t/'index');verified=time.monotonic()
   run=subprocess.run([a.poc,str(t/'index'),str(t/'ref.fa'),str(d/'guides.fa'),'5','2','.75','GG','1'],capture_output=True,text=True,check=True,timeout=90);native=time.monotonic()
   m=re.search(r'states=(\d+) search_wall_seconds=([\deE.+-]+)',run.stderr);assert m
   load=re.search(r'index_load_seconds=([\deE.+-]+)',run.stderr)
   verify(t/'ref.fa',t/'index');reverified=time.monotonic()
   results.append({'region':region['name'],'run':rep,'copy_seconds':copied-start,'verify_before_seconds':verified-copied,'native_total_seconds':native-verified,'native_search_seconds':float(m[2]),'index_constructor_seconds':float(load[1]) if load else None,'native_nonsearch_seconds':native-verified-float(m[2]),'verify_after_seconds':reverified-native,'states':int(m[1])})
   Path(a.out).write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results[-1]),flush=True)
