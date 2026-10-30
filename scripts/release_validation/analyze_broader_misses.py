#!/usr/bin/env python3
"""Join oracle, candidate, imported interval, and direct WFA diagnostic evidence."""
import argparse
import json
import re
from pathlib import Path
from resolve_bulge_discrepancies import fasta, identity, replay, overlaps
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--panels',required=True);p.add_argument('--wfa-log',required=True);p.add_argument('--output',required=True);a=p.parse_args()
raw={}; verified={}
for line in Path(a.wfa_log).read_text().splitlines():
 f=line.split('\t')
 if f[0]=='WFA_RAW': raw[f[1]]={'cigar':f[2],'matches':int(f[3]),'mismatches':int(f[4]),'gaps':int(f[5]),'max_gap_size':int(f[6])}
 if f[0]=='WFA_PROBE':verified[f[1]]=f[2]
records=[]
for panel in ('synthetic','real'):
 d=Path(a.panels)/panel;review=json.loads((d/'alignment_review.json').read_text());probes=json.loads((d/'missing_probes.json').read_text())
 refs, guides = fasta(d/'ref.fa'), fasta(d/'guides.fa')
 exact = []
 for line in (d/'exact.sam').read_text().splitlines():
  if line.startswith('@'):continue
  f=line.split('\t');flag=int(f[1])
  if flag&4:continue
  span=sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]) if op in 'MDN=X')
  start=int(f[3])-1; key=(f[0],f[2],start,start+span,'-' if flag&16 else '+')
  info=identity(key,refs);info.update(cigar_sam=f[5],tags=f[11:],alignment=replay(guides[key[0]],info['oriented_target'],f[5],key[4]=='-'))
  exact.append((key,info))
 for locus in review['metrics']['missing_loci']:
  g,c,strand,ps,pe,pam=locus
  discrepancies=[x for x in review['discrepancies'] if [x['guide'],x['contig'],x['strand'],x['pam_start'],x['pam_end'],x['pam']]==locus]
  pp=[x for x in probes if x['pam_locus']==locus]
  candidates={}
  for row in discrepancies:
   for hit in row['overlapping_raw']:
    candidates[(hit['start'],hit['end'],hit['cigar_sam'])]=hit
  for key,hit in exact:
   if any(overlaps(key,(x['guide'],x['contig'],x['start'],x['end'],x['strand'])) for x in discrepancies):
    candidates[(hit['start'],hit['end'],hit['cigar_sam'])]=hit
  entries=[]
  for v in pp:
   name,contig,start,end,st=v['key']; ident='|'.join(map(str,(panel,name,contig,start,end,st)))
   if ident not in raw or ident not in verified: raise ValueError('missing WFA diagnostic: '+ident)
   entries.append({'interval':v['key'],'wfa_raw':raw[ident],'wfa_accepted':verified[ident]!='None','wfa_result':verified[ident], 'imported_probe':v['attempts']})
  accepted=[e for e in entries if e['wfa_accepted']]
  dist=min((abs(e['interval'][2]-r['start']) for e in accepted for r in candidates.values()),default=None)
  if not accepted: diagnosis='all_oracle_intervals_rejected_by_selected_wfa_alignment'
  elif dist is not None and dist>2: diagnosis='accepted_intervals_start_beyond_two_base_candidate_expansion'
  else: diagnosis='unresolved_candidate_or_expansion'
  records.append({'panel':panel,'locus':locus,'oracle_alignments':discrepancies,
                  'overlapping_raw_candidates':list(candidates.values()),'tested_intervals':entries,
                  'min_start_distance_for_wfa_accepted':dist,'diagnosis':diagnosis})
Path(a.output).write_text(json.dumps(records,indent=2)+'\n')
from collections import Counter
print(json.dumps({'sites':len(records),'diagnoses':dict(Counter(x['diagnosis'] for x in records))}))
