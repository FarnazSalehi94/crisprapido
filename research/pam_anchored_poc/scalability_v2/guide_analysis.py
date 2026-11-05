#!/usr/bin/env python3
"""Characterize oracle-validated panel guides and reported reference-N spans."""
import argparse,collections,json,re,sys
from pathlib import Path
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'research/pam_anchored_poc/validation'));from index_check import fasta
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--guides',required=True);p.add_argument('--paf',required=True);p.add_argument('--out',required=True);a=p.parse_args();refs=dict(fasta(a.reference));guides=dict(fasta(a.guides));by=collections.defaultdict(list)
for line in Path(a.paf).read_text().splitlines():
 f=line.split('\t');lo,hi=int(f[7]),int(f[8]);t=refs[f[5]][lo:hi];by[f[0]].append({'contig':f[5],'strand':f[4],'anchor':hi if f[4]=='+' else lo,'contains_N':'N' in t})
rows=[]
for gid,g in guides.items():
 ops=by[gid];rows.append({'guide':gid,'sequence':g,'gc_fraction':(g.count('G')+g.count('C'))/len(g),'unique_dinucleotides':len({g[i:i+2] for i in range(len(g)-1)}),'max_homopolymer':max(len(m.group()) for m in re.finditer(r'(.)\1*',g)),'sites':len(ops),'N_overlapping_sites':sum(x['contains_N'] for x in ops),'strands':sorted(set(x['strand'] for x in ops))})
result={'guides':rows,'guide_count':len(rows),'gc_min':min(x['gc_fraction'] for x in rows),'gc_max':max(x['gc_fraction'] for x in rows),'unique_dinucleotide_min':min(x['unique_dinucleotides'] for x in rows),'unique_dinucleotide_max':max(x['unique_dinucleotides'] for x in rows),'reported':sum(x['sites'] for x in rows),'N_overlapping_sites':sum(x['N_overlapping_sites'] for x in rows),'N_overlapping_by_strand':dict(collections.Counter(x['strand'] for values in by.values() for x in values if x['contains_N']))}
Path(a.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='guides'}))
