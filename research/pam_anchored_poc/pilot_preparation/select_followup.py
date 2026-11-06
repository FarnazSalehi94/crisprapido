#!/usr/bin/env python3
"""Choose previous repetitive cases and bounded windows around real N runs."""
import argparse, json, mmap, re, sys
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--stage',required=True); p.add_argument('--fai',required=True); p.add_argument('--out',required=True); p.add_argument('--max-runs',type=int,default=3000); a=p.parse_args()
source=Path(__file__).resolve().parents[1]/'large_coordinate/select_probes.py'
# Reuse only input validation and window extraction, without its selection or guide search.
sys.argv=[str(source),'--stage',a.stage,'--fai',a.fai,'--out',a.out]
ns={};exec(source.read_text().split('cases=[]')[0],ns)
rows=ns['rows'];positions=ns['positions'];candidates=ns['candidates'];cases=[]
previous=json.loads((source.parent/'results/cases.json').read_text())['cases'];skip=json.loads((source.parent/'results/main-summary.json').read_text())['skips'];follow=json.loads((source.parent/'results/followup-case.json').read_text())['cases'][0]
for x in skip:
 c=next(c for c in previous if all(c[k]==x[k] for k in ('label','strand','global_start')))
 if c['global_start']!=follow['global_start'] or c['strand']!=follow['strand']:cases.append(c)
assert len(cases)==13
wanted={(kind,strand) for kind in ('target_n','pam_required_n','pam_wildcard_n') for strand in '+-'}; boundary_wanted=set('+-');runs=0
for i,r in enumerate(rows):
 with (Path(a.stage)/'reference.fa').open('rb') as fh:
  with mmap.mmap(fh.fileno(),0,access=mmap.ACCESS_READ) as mm:
   cursor=r['byte_offset'];stop=(rows[i+1]['byte_offset'] if i+1<len(rows) else len(mm))
   while cursor<stop and runs<a.max_runs and (wanted or boundary_wanted):
    pos=mm.find(b'N',cursor,stop)
    if pos<0:break
    # Include both ends; FASTA line breaks may split an N run.
    match=re.compile(b'[ACGT]').search(mm,pos,stop); end=match.start() if match else stop; runs+=1
    for byte in (pos,max(pos,end-1)):
     delta=byte-r['byte_offset']; center=(delta//r['line_width'])*r['line_bases']+min(delta%r['line_width'],r['line_bases']-1)
     lo=max(0,center-80);hi=min(r['length'],center+80)
     if hi-lo<23:continue
     for kind,strand in sorted(wanted | {('target_n',s) for s in boundary_wanted if center<1000 or center>r['length']-1000}):
      possible=[x for x in candidates(i,lo,hi,strand,kind) if x[3][:20].count('N')<=5]
      if not possible:continue
      start,g,t,orient,shadow=possible[0];anchor=start+20 if strand=='+' else start+3
      label=('boundary_target_n' if kind=='target_n' and (center<1000 or center>r['length']-1000) else 'natural_'+kind)
      case={'label':label,'kind':kind,'contig':r['name'],'contig_index':i,'contig_length':r['length'],'global_start':positions[i]+start,'start':start,'end':start+23,'strand':strand,'anchor':anchor,'reference_genomic':g,'indexed_genomic':t,'oriented_reference':orient,'oriented_index':shadow,'guide':shadow[:20],'pam':orient[20:],'true_mismatches':sum(x!=y for x,y in zip(shadow[:20],orient[:20])),'strict_pam':orient[20] in 'ACGT' and orient[21:]=='GG'}
      if (kind,strand) in wanted:cases.append(case);wanted.remove((kind,strand))
      if label=='boundary_target_n' and strand in boundary_wanted:
       if case not in cases:cases.append(case)
       boundary_wanted.remove(strand)
    cursor=end+1
   if runs>=a.max_runs or (not wanted and not boundary_wanted):break
Path(a.out).write_text(json.dumps({'cases':cases,'n_runs_inspected':runs,'unfound':sorted(wanted),'boundary_strands_unfound':sorted(boundary_wanted),'selection':'bounded scan of N runs, exact pattern probes only'},indent=2)+'\n')
print('cases',len(cases),'N runs',runs,'unfound',sorted(wanted),'boundary',sorted(boundary_wanted))
