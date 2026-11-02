#!/usr/bin/env python3
"""Audit outputs and summarize per-guide scaling without altering search results."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
p=argparse.ArgumentParser(description=__doc__)
for name in ('panel','results','samtools','matrices','output'):p.add_argument('--'+name,required=True)
a=p.parse_args();panel=Path(a.panel).resolve();results=Path(a.results).resolve();rows=json.loads((results/'results.json').read_text());out=Path(a.output);out.mkdir(exist_ok=False)
script=Path(__file__).resolve().parent/'audit_emitted.py';prepared=set();summary=[]
for r in rows:
 if r['status']!='PASS':raise ValueError('nonpassing search')
 win=panel/r['window'];run=results/r['window']/r['guide'];ref=win/'ref.fa'
 if r['window'] not in prepared:
  subprocess.run([a.samtools,'faidx',str(ref)],check=True);prepared.add(r['window'])
 audit=out/(r['window']+'__'+r['guide']);subprocess.run([sys.executable,str(script),'--paf',str(run/'actual.paf'),'--fasta',str(ref),'--guides',str(run/'guide.fa'),'--matrices',a.matrices,'--samtools',a.samtools,'--output',str(audit)],stdout=subprocess.DEVNULL,check=True)
 text=(run/'time.txt').read_text();line=next(x for x in text.splitlines() if 'Elapsed (wall clock)' in x);v=line.rsplit(': ',1)[1].split(':');elapsed=sum(float(x)*60**i for i,x in enumerate(reversed(v)))
 profile=r['profile'];summary.append({'window':r['window'],'guide':r['guide'],'time_seconds':elapsed,'peak_rss_kib':r['peak_rss_kib'],
 'raw_e5':r['raw_e5'],'raw_e0':r['raw_e0'],'unique_raw':r['unique_raw'],'wfa_calls':r['wfa_calls'],
 'wfa_per_candidate':r['wfa_per_unique_candidate'],'reported':r['reported_rows'],'audit':'PASS',
 'columba_seconds':profile['columba_wait']['nanoseconds']/1e9,'expansion_seconds':profile['expansion_enumeration']['nanoseconds']/1e9,
 'wfa_seconds':profile['wfa_align']['nanoseconds']/1e9,'alternative_seconds':profile['eligible_enumeration']['nanoseconds']/1e9,
 'pam_seconds':profile['pam_extract']['nanoseconds']/1e9,'cfd_seconds':profile['cfd_score']['nanoseconds']/1e9,'selection_seconds':profile['representative_selection']['nanoseconds']/1e9})
with (out/'per_guide.tsv').open('x',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
