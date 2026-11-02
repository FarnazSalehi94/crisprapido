#!/usr/bin/env python3
"""Count candidate keys and same-reference-span verification intervals."""
import argparse,json,re
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sam',action='append',required=True);p.add_argument('--output',required=True);a=p.parse_args();result=[]
for name in a.sam:
 path=Path(name);intervals=set();keys=set();raw=0
 for line in path.open():
  if line.startswith('@'):continue
  f=line.split('\t');flag=int(f[1])
  if flag&4:continue
  pos=int(f[3])-1;span=sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]) if op in 'MDN=X');raw+=1
  keys.add((f[0],f[2],pos,bool(flag&16),f[5]));intervals.add((f[0],f[2],pos,pos+span,bool(flag&16)))
 result.append({'sam':str(path),'raw':raw,'candidate_keys':len(keys),'verification_intervals':len(intervals),'repeated_verification_intervals':len(keys)-len(intervals)})
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
