#!/usr/bin/env python3
"""Compare PAM-anchor sets and representative alignments, not only CIGAR identity."""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--baseline',required=True);p.add_argument('--candidate',required=True);p.add_argument('--out',required=True);a=p.parse_args()
def read(path):
 rows={};duplicates=[]
 for line in Path(path).read_text().splitlines():
  f=line.split('\t');assert len(f)>=12 and f[0]=='pilot' and f[4] in ('+','-');tags={z.split(':',2)[0]:z.split(':',2)[2] for z in f[12:]};anchor=int(f[8]) if f[4]=='+' else int(f[7])-3;key=(f[0],f[5],f[4],anchor)
  if key in rows:duplicates.append(key)
  rows[key]={'start':int(f[7]),'end':int(f[8]),'cigar':tags['cg'],'mismatches':tags['nm'],'events':tags['ng'],'size':tags['bs'],'cfd':tags.get('cf')}
 return rows,duplicates
old,olddup=read(a.baseline);new,newdup=read(a.candidate);common=old.keys()&new.keys();representative=[k for k in sorted(common) if {x:y for x,y in old[k].items() if x!='cfd'}!={x:y for x,y in new[k].items() if x!='cfd'}];scores=[k for k in sorted(common) if old[k]['cfd']!=new[k]['cfd']];known=('pilot','HG00735#1#JAHBCH020000014.1','+',43957024)
result={'baseline_sites':len(old),'candidate_sites':len(new),'shared':len(common),'missing':len(old.keys()-new.keys()),'unexpected':len(new.keys()-old.keys()),'baseline_duplicate_anchors':len(olddup),'candidate_duplicate_anchors':len(newdup),'representative_alignment_differences':len(representative),'cfd_differences':len(scores),'known_exact_target_recovered':known in new and new[known]['cigar']=='20=' and new[known]['mismatches']=='0','alignment_difference_examples':[{'key':k,'baseline':old[k],'candidate':new[k]} for k in representative[:10]],'cfd_difference_examples':[{'key':k,'baseline':old[k]['cfd'],'candidate':new[k]['cfd']} for k in scores[:10]]}
Path(a.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));
if olddup or newdup or result['missing'] or result['unexpected'] or len(old)!=21447 or not result['known_exact_target_recovered']:raise SystemExit(1)
