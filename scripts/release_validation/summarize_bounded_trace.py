#!/usr/bin/env python3
"""Tabulate per-site candidate geometry, oracle alignment and WFA2 selection."""
import argparse
import csv
import json
import re
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--prior',required=True);p.add_argument('--corrected',required=True);p.add_argument('--output',required=True);a=p.parse_args()
prior=json.loads(Path(a.prior).read_text()); rows=[]
for item in prior:
 panel=item['panel']; guide,contig,strand,ps,pe,pam=item['locus']
 accepted=[x for x in item['tested_intervals'] if x['wfa_accepted']]
 choices=accepted or item['tested_intervals']
 pairs=[(abs(e['interval'][2]-c['start']),e,c) for e in choices for c in item['overlapping_raw_candidates']]
 _,entry,candidate=min(pairs,key=lambda x:(x[0],x[1]['interval'][2],x[2]['start']))
 s,e=entry['interval'][2:4]
 oracle=[v for x in item['oracle_alignments'] if x['start']==s and x['end']==e for v in x['oracle_alignments']]
 if not oracle: raise ValueError('missing oracle alignment')
 def cost(v):return 3*v['mismatches'] + (6 if v['bulge']!='none' else 0)
 chosen=min(oracle,key=cost)
 new=[]
 for line in (Path(a.corrected)/panel/'out.paf').read_text().splitlines():
  f=line.split('\t');anchor=int(f[8]) if strand=='+' else int(f[7]); desired=ps if strand=='+' else pe
  if f[0]==guide and f[5]==contig and f[4]==strand and anchor==desired:
   new.append(f)
 if len(new)!=1:raise ValueError('expected exactly one corrected PAM anchor')
 w=entry['wfa_raw'];final=new[0]
 rows.append(dict(panel=panel,guide=guide,contig=contig,strand=strand,pam_start=ps,pam_end=pe,pam=pam,
  raw_start=candidate['start'],raw_end=candidate['end'],raw_sam_cigar=candidate['cigar_sam'],
  target_start=s,target_end=e,start_offset=s-candidate['start'],
  oracle_cigar=chosen['cigar'],oracle_mismatches=chosen['mismatches'],
  oracle_edit_distance=chosen['mismatches'] + (chosen['bulge']!='none'), oracle_score_cost=cost(chosen),
  selected_wfa_cigar=w['cigar'],selected_wfa_mismatches=w['mismatches'],
  selected_wfa_edit_distance=w['cigar'].count('X')+w['cigar'].count('I')+w['cigar'].count('D'),
  selected_wfa_score_cost=3*w['cigar'].count('X')+sum(5+len(x) for x in re.findall(r'I+|D+',w['cigar'])),
  selected_wfa_gaps=w['gaps'],selected_wfa_max_gap_size=w['max_gap_size'],
  old_diagnosis=item['diagnosis'],corrected_start=final[7],corrected_end=final[8],
  corrected_cigar=next(x[5:] for x in final[12:] if x.startswith('cg:Z:'))))
with open(a.output,'w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
