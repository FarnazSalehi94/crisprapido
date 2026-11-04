#!/usr/bin/env python3
"""Independently replay and score both representatives at every shared anchor."""
import csv,json,re,subprocess,sys
from collections import Counter
from pathlib import Path
repo=Path(__file__).resolve().parents[3];sys.path.insert(0,str(repo/'scripts/release_validation'));from audit_emitted import independent_cfd,matrices
sys.path.insert(0,str(repo/'research/pam_anchored_poc/validation'));from index_check import fasta
ref,guide,poc,baseline,probe,out=sys.argv[1:];refs=dict(fasta(ref));guides=dict(fasta(guide));mm,pm=matrices(repo);rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def replay(gid,c,strand,start,end,cigar):
 seq=refs[c];g=guides[gid];assert 0<=start<end<=len(seq);t=seq[start:end] if strand=='+' else rc(seq[start:end]);pam=seq[end:end+3] if strand=='+' else rc(seq[start-3:start]);assert len(pam)==3 and pam[0] in 'ACGT' and pam[1:]=='GG'
 ops=[(int(n),op) for n,op in re.findall(r'(\d+)([=XID])',cigar)];assert ''.join(str(n)+op for n,op in ops)==cigar
 qi=ti=matches=mis=events=gaplen=maxgap=0;prev='';a=[];b=[];gaps=[]
 for n,op in ops:
  if op in 'ID':
   if not prev or prev not in 'ID':events+=1;gaplen=0
   gaplen+=n;maxgap=max(maxgap,gaplen);gaps.append({'operation':op,'guide_offset':qi,'target_offset':ti,'length':n})
  else:gaplen=0
  prev=op
  for _ in range(n):
   if op=='I':a.append(g[qi]);b.append('-');qi+=1
   elif op=='D':a.append('-');b.append(t[ti]);ti+=1
   else:
    assert qi<len(g) and ti<len(t);same=g[qi]==t[ti];assert same==(op=='=');a.append(g[qi]);b.append(t[ti]);qi+=1;ti+=1;matches+=same;mis+=not same
 assert qi==len(g) and ti==len(t) and mis<=5 and events<=1 and maxgap<=2 and matches>=15
 return {'start':start,'end':end,'cigar':cigar,'pam':pam,'mismatches':mis,'gap_events':events,'max_gap_size':maxgap,'gap_placements':gaps,'matches':matches,'match_fraction':matches/len(g),'affine_penalty':3*mis+sum(5+n for n,op in ops if op in 'ID'),'cfd':independent_cfd(''.join(a),''.join(b),ops,pam,mm,pm),'_target':t}
def key(g,c,s,a):return (g,c,s,a)
base={}
for line in Path(baseline).read_text().splitlines():
 f=line.split('\t');g,c,s=f[0],f[5],f[4];start,end=int(f[7]),int(f[8]);tags={x.split(':',2)[0]:x.split(':',2)[2] for x in f[12:]};k=key(g,c,s,end if s=='+' else start);assert k not in base;base[k]=replay(g,c,s,start,end,tags['cg']);base[k]['reported_cfd']=float(tags['cf']) if 'cf' in tags else None
 if base[k]['reported_cfd'] is not None:assert base[k]['cfd'] is not None and abs(base[k]['reported_cfd']-base[k]['cfd'])<=.000051
lines=Path(poc).read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'));rows=list(csv.DictReader(lines[i:],delimiter='\t'));new={}
for r in rows:
 k=key(r['guide'],r['contig'],r['strand'],int(r['anchor']));assert k not in new;new[k]=replay(r['guide'],r['contig'],r['strand'],int(r['start']),int(r['end']),r['cigar']);assert new[k]['mismatches']==int(r['mismatches']) and new[k]['gap_events']==int(r['bulge_events']) and new[k]['matches']==int(r['match_count'])
assert set(base)==set(new)
differences=[];queries=[];scores=[]
for k in sorted(base):
 a=base[k];b=new[k]
 if (a['start'],a['end'],a['cigar'])==(b['start'],b['end'],b['cigar']):continue
 for r in (a,b):queries.append('\t'.join((guides[k[0]],r['_target'],r['cigar'],r['pam'][1:])));scores.append(r['cfd'])
 a=dict(a);b=dict(b);a.pop('_target');b.pop('_target');different_cfd=(a['cfd'] is None)!=(b['cfd'] is None) or (a['cfd'] is not None and b['cfd'] is not None and abs(a['cfd']-b['cfd'])>1e-9)
 differences.append({'key':list(k),'baseline':a,'prototype':b,'mismatch_count_differs':a['mismatches']!=b['mismatches'],'bulge_event_count_differs':a['gap_events']!=b['gap_events'],'bulge_max_size_differs':a['max_gap_size']!=b['max_gap_size'],'coordinates_differ':(a['start'],a['end'])!=(b['start'],b['end']),'bulge_placements_differ':a['gap_placements']!=b['gap_placements'],'cfd_differs':different_cfd,'affine_preference':('baseline' if a['affine_penalty']<b['affine_penalty'] else 'prototype' if b['affine_penalty']<a['affine_penalty'] else 'tie'), 'preference':'No biological preference established; select under a declared deterministic alignment objective.'})
r=subprocess.run([probe],input='\n'.join(queries)+'\n',capture_output=True,text=True,check=True,cwd=repo);actual=r.stdout.splitlines();assert len(actual)==len(scores)
for a,b in zip(actual,scores):assert (a=='NA' and b is None) or (a!='NA' and b is not None and abs(float(a)-b)<1e-9),(a,b)
summary={'total_anchors':len(base),'differences':len(differences),'both_independently_eligible':len(differences),'mismatch_count_differs':sum(d['mismatch_count_differs'] for d in differences),'bulge_event_count_differs':sum(d['bulge_event_count_differs'] for d in differences),'bulge_max_size_differs':sum(d['bulge_max_size_differs'] for d in differences),'coordinates_differ':sum(d['coordinates_differ'] for d in differences),'bulge_placements_differ':sum(d['bulge_placements_differ'] for d in differences),'cfd_differs':sum(d['cfd_differs'] for d in differences),'affine_preference_counts':dict(Counter(d['affine_preference'] for d in differences)),'rust_independent_scores_compared':len(scores)}
Path(out).mkdir(exist_ok=True,parents=True);Path(out,'representatives.json').write_text(json.dumps(differences,indent=2)+'\n');Path(out,'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,sort_keys=True))
