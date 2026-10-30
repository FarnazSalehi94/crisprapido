#!/usr/bin/env python3
"""Summarize the previously missing PAM anchors and recovered final alignments."""
import argparse
import csv
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline', required=True)
p.add_argument('--corrected', required=True)
p.add_argument('--output', required=True)
a=p.parse_args()
rows=[]
for panel in ('indexed_synthetic','indexed_real_chm13','indexed_12_windows'):
    old=json.loads((Path(a.baseline)/panel/'alignment_review.json').read_text())
    missing={tuple(x) for x in old['metrics']['missing_loci']}
    recovered={}
    for line in (Path(a.corrected)/panel/'out.paf').read_text().splitlines():
        f=line.split('\t'); anchor=int(f[8]) if f[4]=='+' else int(f[7]);
        bounds=(anchor,anchor+3) if f[4]=='+' else (anchor-3,anchor)
        key=(f[0],f[5],f[4],bounds[0],bounds[1])
        recovered[key]=(f[7],f[8],next(x[5:] for x in f[12:] if x.startswith('cg:Z:')))
    for locus in sorted(missing):
        candidates=[x for x in old['discrepancies'] if x['classification']=='C' and
                    (x['guide'],x['contig'],x['strand'],x['pam_start'],x['pam_end'],x['pam'])==locus]
        if not candidates: raise ValueError('missing oracle evidence')
        first=candidates[0]; raw=first['overlapping_raw']
        if len(raw)!=1: raise ValueError('ambiguous original raw candidate')
        rk=(locus[0],locus[1],locus[2],locus[3],locus[4])
        if rk not in recovered: raise ValueError('corrected locus not reported')
        start,end,cigar=recovered[rk]
        rows.append(dict(panel=panel,guide=locus[0],contig=locus[1],strand=locus[2],
            oracle_start=first['start'],oracle_end=first['end'],pam_start=locus[3],pam_end=locus[4],pam=locus[5],
            original_raw_start=raw[0]['start'],original_raw_end=raw[0]['end'],original_raw_cigar=raw[0]['cigar_sam'],
            original_raw_alignment=raw[0]['alignment']['normalized_cigar'],original_raw_pam=raw[0]['pam'],
            original_expansion='not_triggered_on_accepted_primary', original_pam_in_final='absent',
            corrected_start=start,corrected_end=end,corrected_cigar=cigar))
with open(a.output,'w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
