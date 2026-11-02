#!/usr/bin/env python3
"""Compare an earlier small SAM's keys to a larger saved SAM, streaming the latter."""
import argparse,json,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('old','new','output'):p.add_argument('--'+n,required=True)
a=p.parse_args();t=time.monotonic()
def key(line):
 f=line.split('\t',6);flag=int(f[1])
 return None if flag&4 else (f[0],f[2],int(f[3])-1,bool(flag&16),f[5])
old=set()
with open(a.old) as f:
 for line in f:
  if not line.startswith('@'):
   k=key(line)
   if k:old.add(k)
oldstarts={(q,c,s,r) for q,c,s,r,cigar in old};found=set();foundstarts=set();n=0
with open(a.new) as f:
 for line in f:
  if line.startswith('@'):continue
  k=key(line)
  if not k:continue
  n+=1
  if k in old:found.add(k)
  start=k[:4]
  if start in oldstarts:foundstarts.add(start)
missing=old-found;missingstart=oldstarts-foundstarts
result={'old_unique_keys':len(old),'new_mapped_rows':n,'exact_key_intersection':len(found),'old_keys_not_exactly_present':len(missing),'old_distinct_starts':len(oldstarts),'old_starts_present_in_new':len(foundstarts),'old_starts_not_present_in_new':len(missingstart),'old_keys_absent_examples':[list(x) for x in sorted(missing)[:20]],'old_starts_absent_examples':[list(x) for x in sorted(missingstart)[:20]],'elapsed_seconds':time.monotonic()-t,'limitation':'A different starting coordinate may overlap the same eligible PAM locus; this measures candidate keys and starts, not recall.'}
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
