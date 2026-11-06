#!/usr/bin/env python3
"""Validate bounded probe logs against selected cases and independent source windows."""
import argparse,collections,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--cases',required=True);p.add_argument('--independent',required=True);p.add_argument('--run',required=True);p.add_argument('--cap',type=int,required=True);p.add_argument('--out',required=True);a=p.parse_args();cases=json.loads(Path(a.cases).read_text())['cases'];ind=json.loads(Path(a.independent).read_text())['cases'];run=Path(a.run);measure=json.loads((run/'measurements.json').read_text());assert measure['returncode']==0 and not measure['stop_reason'] and 'PHASE\tpostverification' in (run/'stderr.log').read_text() and 'PHASE\tdone' in (run/'stderr.log').read_text()
expected={(c['label'],c['strand'],c['global_start']):c for c in cases};assert len(expected)==len(cases);independent={(c['label'],c['strand'],c['global_start']):c for c in ind};assert set(expected)==set(independent)
seen={};controls=[]
for line in (run/'stdout.log').read_text().splitlines():
 f=line.split('\t')
 if f[0]=='CASE':
  assert len(f)==11;k=(f[1],f[3],int(f[4]));assert k in expected and k not in seen;c=expected[k];assert f[2]==c['kind'] and f[5]==c['contig'] and int(f[6])==c['start'] and int(f[7])==c['anchor'] and 1<=int(f[8])<=a.cap
  assert (int(f[9]),int(f[10]))==((0,-1) if c['kind']=='cross_boundary' else (1,c['true_mismatches']));seen[k]={'result':'pass','range_width':int(f[8])}
 elif f[0]=='SKIP':
  assert len(f)==6;k=(f[1],f[3],int(f[4]));assert k in expected and k not in seen and int(f[5])>a.cap;seen[k]={'result':'skipped_cap','range_width':int(f[5])}
 elif f[0]=='CONTROL':
  assert len(f)==6;assert f[1] in ('target_n','pam_wildcard_n','pam_required_n') and f[2] in '+-';assert (int(f[4]),int(f[5]))==((1,1) if f[1]=='target_n' else (0,-1));controls.append({'kind':f[1],'strand':f[2],'global_start':int(f[3]),'emitted':int(f[4]),'mismatches':int(f[5])})
assert set(seen)==set(expected)
result={'selected':len(cases),'completed':sum(x['result']=='pass' for x in seen.values()),'skipped':sum(x['result']=='skipped_cap' for x in seen.values()),'completed_by_label_strand':{},'skips':[],'controls':controls}
for k,v in seen.items():
 if v['result']=='pass':key=k[0]+' '+k[1];result['completed_by_label_strand'][key]=result['completed_by_label_strand'].get(key,0)+1
 else:result['skips'].append({'label':k[0],'strand':k[1],'global_start':k[2],'range_width':v['range_width']})
Path(a.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='controls'}))
