#!/usr/bin/env python3
"""Sample the job shell's Linux process tree; summed RSS may double-count pages."""
import os,sys,time
from pathlib import Path
root=int(sys.argv[1]);peak=0
while True:
 records={}
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   lines=(p/'status').read_text().splitlines();d={x.split(':',1)[0]:x.split(':',1)[1].strip() for x in lines if ':' in x};records[int(p.name)]=(int(d['PPid']),int(d.get('VmRSS','0 kB').split()[0]))
  except (OSError,KeyError,ValueError):pass
 if root not in records:break
 tree={root}
 while True:
  expanded=tree|{pid for pid,(parent,rss) in records.items() if parent in tree}
  if expanded==tree:break
  tree=expanded
 rss=sum(records[pid][1] for pid in tree if pid in records);peak=max(peak,rss)
 print('%s members=%d summed_rss_kib=%d sampled_peak_kib=%d'%(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),len(tree),rss,peak),flush=True);time.sleep(1)
