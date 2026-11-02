#!/usr/bin/env python3
"""Check whether saved raw candidates overlap prior reported reference spans."""
import argparse,bisect,collections,json,re,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('old-paf','new-sam','output'):p.add_argument('--'+n,required=True)
a=p.parse_args();t=time.monotonic();regions=collections.defaultdict(list)
with open(a.old_paf) as stream:
 for line in stream:
  f=line.split('\t');key=(f[0],f[5],f[4]);start,end=int(f[7]),int(f[8])
  if not 0 < end-start <= 21:raise ValueError('prior interval outside expected bound')
  regions[key].append((start,end,line.rstrip('\n')))
starts={k:[x[0] for x in sorted(v)] for k,v in regions.items()};regions={k:sorted(v) for k,v in regions.items()};hits=set();n=0
with open(a.new_sam) as stream:
 for line in stream:
  if line.startswith('@'):continue
  f=line.split('\t',6);flag=int(f[1])
  if flag&4:continue
  n+=1;pos=int(f[3])-1;k=(f[0],f[2],'-' if flag&16 else '+')
  if k not in regions:continue
  # A 20-base query and edit bound 7 imply reference span <=27 bases.
  # The old reported b=1/z=1 intervals have span <=21 bases.
  span=sum(int(v) for v,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]) if op in 'MDN=X')
  if span>27:raise ValueError('candidate exceeds assumed 27-base bound')
  lo=bisect.bisect_left(starts[k],pos-20);hi=bisect.bisect_left(starts[k],pos+span)
  if lo==hi:continue
  for oldstart,oldend,row in regions[k][lo:hi]:
   if oldstart<pos+span and oldend>pos:hits.add(row)
rows=[x[2] for v in regions.values() for x in v];missing=[r for r in rows if r not in hits]
result={'old_reported_rows':len(rows),'new_sam_mapped_rows':n,'old_rows_with_overlapping_candidate':len(hits),'old_rows_without_overlapping_candidate':len(missing),'missing':missing,'elapsed_seconds':time.monotonic()-t,'limitation':'Same-guide/contig/strand interval overlap is a necessary candidate relation only; it neither proves equivalence nor guarantees final recovery.'}
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
