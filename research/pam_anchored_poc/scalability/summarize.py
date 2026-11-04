#!/usr/bin/env python3
"""Replay each run, check repeat determinism and optionally enumerate 100 kb oracle."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--dataset',required=True);p.add_argument('--runs',required=True);p.add_argument('--out',required=True);a=p.parse_args();data=Path(a.dataset);runs=Path(a.runs);out=Path(a.out);out.mkdir(exist_ok=True,parents=True)
regions=json.loads((data/'regions.json').read_text());measurements=json.loads((runs/'runs.json').read_text());result=[]
for region in regions:
 name=region['name'];rs=[x for x in measurements if x['region']==name];d=data/name;r=runs/name
 hashes={}; successful={}
 for backend in ('columba-wfa2','pam-anchored'):
  files=[r/(backend+'-'+str(x['run'])+'.paf') for x in rs if x['backend']==backend and x['status']==0]
  hashes[backend]=[hashlib.sha256(f.read_bytes()).hexdigest() for f in files];successful[backend]=files
 item={'region':name,'runs':rs,'output_hashes':hashes,'deterministic':{k:len(set(v))<=1 for k,v in hashes.items()}}
 if all(successful.values()):
  path=out/(name+'.json');cmd=[sys.executable,str(Path(__file__).with_name('audit.py')),'--reference',str(d/'ref.fa'),'--guides',str(d/'guides.fa'),'--baseline',str(successful['columba-wfa2'][0]),'--experimental',str(successful['pam-anchored'][0]),'--out',str(path)]
  if region['bases']==100000:cmd.append('--oracle')
  subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL,timeout=180)
  audit=json.loads(path.read_text());item['audit']={k:v for k,v in audit.items() if k!='alignments'}
  for k,v in hashes.items():
   if len(set(v))>1:raise RuntimeError('nondeterministic output; audit all runs: '+name+' '+k)
 result.append(item);(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(name,json.dumps(item.get('audit',{})),flush=True)
