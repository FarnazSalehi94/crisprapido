#!/usr/bin/env python3
"""Simple repeat and nearby-PAM proxies; no external annotation claimed."""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--dataset',required=True);p.add_argument('--out',required=True);a=p.parse_args();root=Path(a.dataset);results=[]
for meta in json.loads((root/'regions.json').read_text()):
 lines=(root/meta['name']/'ref.fa').read_text().splitlines();seq=''.join(lines[1:]).upper();positions=[(len(seq)-20)*i//99 for i in range(100)];sampled=repeated=0
 for i in positions:
  k=seq[i:i+20]
  if set(k)-set('ACGT'):continue
  sampled+=1
  if (i>=20 and seq.find(k,0,i-19)!=-1) or seq.find(k,i+20)!=-1:repeated+=1
 counts={};near={}
 for strand in ('+','-'):
  positions=[i for i in range(len(seq)-2) if (seq[i] in 'ACGT' and seq[i+1:i+3]=='GG') if strand=='+'] if strand=='+' else [i+3 for i in range(len(seq)-2) if seq[i:i+2]=='CC' and seq[i+2] in 'ACGT']
  counts[strand]=len(positions);near[strand]=sum(b-a<=3 for a,b in zip(positions,positions[1:]))
 results.append({'region':meta['name'],'sampled_resolved_20mers':sampled,'repeated_nonoverlap_20mers':repeated,'pam_counts':counts,'adjacent_pam_pairs_at_most_3bp':near})
Path(a.out).write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
