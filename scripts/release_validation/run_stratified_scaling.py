#!/usr/bin/env python3
"""Run separate bounded per-guide regional searches and verify profiler output."""
import argparse
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import time
from multi_exact_scanner import read_guides
p=argparse.ArgumentParser(description=__doc__)
for n in ('panel','output','binary','profile-binary','columba','time-bin'):p.add_argument('--'+n,required=True)
p.add_argument('--timeout',type=int,default=90);p.add_argument('--memory-gib',type=int,default=4)
a=p.parse_args();root=Path(a.panel).resolve();out=Path(a.output).resolve();out.mkdir(exist_ok=False)
windows=json.loads((root/'windows.json').read_text());guides=read_guides(root/'guides.fa');all_results=[]
if a.timeout<1 or a.memory_gib<1:p.error('limits must be positive')
def limit():
 n=a.memory_gib*1024**3;resource.setrlimit(resource.RLIMIT_AS,(n,n))
for win in windows:
 for name,guide in guides.items():
  d=out/win['name']/name;d.mkdir(parents=True);(d/'guide.fa').write_text('>'+name+'\n'+guide+'\n');tmp=d/'tmp';tmp.mkdir()
  window=root/win['name'];common=['-r',str(window/'ref.fa'),'-G',str(d/'guide.fa'),'-p','GG','-m','4','-b','1','-z','1','-f','0.75','-t','1','--columba-bin',a.columba,'--columba-index',str(window/'index')]
  env=dict(os.environ,TMPDIR=str(tmp));stats={'window':win['name'],'guide':name,'region':win['region']}
  try:
   start=time.monotonic()
   with (d/'actual.paf').open('x') as stdout,(d/'actual.stderr').open('x') as stderr:
    subprocess.run([a.time_bin,'-v','-o',str(d/'time.txt'),a.binary]+common+['--keep-columba-sam'],stdout=stdout,stderr=stderr,env=env,preexec_fn=limit,check=True,timeout=a.timeout)
   stats['wall_seconds']=time.monotonic()-start
   start=time.monotonic()
   with (d/'profile.paf').open('x') as stdout,(d/'profile.stderr').open('x') as stderr:
    subprocess.run([a.profile_binary]+common,stdout=stdout,stderr=stderr,env=env,preexec_fn=limit,check=True,timeout=a.timeout)
   stats['profile_wall_seconds']=time.monotonic()-start
   if (d/'actual.paf').read_bytes()!=(d/'profile.paf').read_bytes():raise ValueError('profile output differs')
   profile={}
   for line in (d/'profile.stderr').read_text().splitlines():
    f=line.split('\t')
    if f[0]=='PROFILE':profile[f[1]]={'calls':int(f[2]),'nanoseconds':int(f[3])}
   if 'wfa_align' not in profile:raise ValueError('missing WFA profiler data')
   raw=[];by_e={}
   for sam in tmp.glob('*/output.sam'):
    text=sam.read_text().splitlines();pg=next(x for x in text if x.startswith('@PG\t'))
    m=re.search(r'(?:^| )-e (\d+)(?: |$)',pg)
    if not m:raise ValueError('missing SAM edit bound')
    e=m.group(1); count=0
    for line in text:
     if line.startswith('@'):continue
     f=line.split('\t');flag=int(f[1])
     if flag&4:continue
     count+=1;raw.append((f[0],f[2],f[3],bool(flag&16),f[5]))
    by_e[e]=count
   if set(by_e)!={'0','5'}:raise ValueError('missing retained Columba SAM')
   time_text=(d/'time.txt').read_text();m=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',time_text)
   if not m:raise ValueError('missing RSS')
   stats.update(status='PASS',raw_e5=by_e['5'],raw_e0=by_e['0'],unique_raw=len(set(raw)),
       peak_rss_kib=int(m.group(1)),wfa_calls=profile['wfa_align']['calls'],
       wfa_per_unique_candidate=profile['wfa_align']['calls']/len(set(raw)) if raw else None,
       reported_rows=len((d/'actual.paf').read_text().splitlines()),profile=profile)
  except (subprocess.CalledProcessError,subprocess.TimeoutExpired,ValueError,OSError) as exc:
   stats.update(status='FAIL',error=repr(exc))
  (d/'result.json').write_text(json.dumps(stats,indent=2)+'\n');all_results.append(stats)
  print(win['name'],name,stats['status'],flush=True)
  if stats['status']!='PASS':
   (out/'results.json').write_text(json.dumps(all_results,indent=2)+'\n');raise SystemExit(1)
(out/'results.json').write_text(json.dumps(all_results,indent=2)+'\n')
