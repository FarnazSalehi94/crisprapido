#!/usr/bin/env python3
"""Compare current batch against single-guide calls on a bounded indexed panel."""
import argparse,json,re,subprocess,tempfile,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--cli',required=True);p.add_argument('--native',required=True);p.add_argument('--index',required=True);p.add_argument('--ref',required=True);p.add_argument('--guides',required=True);p.add_argument('--out',required=True);p.add_argument('--repeats',type=int,default=3);a=p.parse_args();guides=[]
for line in Path(a.guides).read_text().splitlines():
 if line.startswith('>'):guides.append([line, ''])
 else:guides[-1][1]+=line.strip()
assert len(guides)==64
base=[a.cli,'--search-backend','pam-anchored','--pam-search-bin',a.native,'--columba-index',a.index,'-r',a.ref,'-t','1','-m','5','-b','1','-z','2','-f','.75','-p','GG'];results=[]
with tempfile.TemporaryDirectory(prefix='crispr-batch-profile-') as td:
 paths=[]
 for i,(name,seq) in enumerate(guides):
  f=Path(td)/('%03d.fa'%i);f.write_text(name+'\n'+seq+'\n');paths.append(f)
 def run(paths):
  start=time.monotonic();lines=[];states=0
  for path in paths:
   r=subprocess.run(base+['-G',str(path)],capture_output=True,text=True,check=True,timeout=60);lines.extend(r.stdout.splitlines());m=re.search(r'states=(\d+)',r.stderr);assert m;states+=int(m.group(1))
  return {'seconds':time.monotonic()-start,'sites':len(lines),'states':states},sorted(lines)
 for i in range(a.repeats):
  batch,b=run([a.guides]);single,s=run(paths);assert b==s and batch['states']==single['states'];results.append({'repetition':i,'batch':batch,'separate':single,'identical_paf_multiset':True})
  Path(a.out).write_text(json.dumps({'guides':len(guides),'runs':results,'scope':'existing 8kb panel; same machine and prebuilt index; does not quantify whole-genome batch memory'},indent=2)+'\n')
print(json.dumps(results))
