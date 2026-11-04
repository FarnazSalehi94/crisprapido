#!/usr/bin/env python3
"""Check all saved difficult panel anchor keys and existing fallback examples."""
import csv, json, sys
from pathlib import Path
root=Path(sys.argv[1]); repo=Path(__file__).resolve().parents[2]
for name in ('indexed_synthetic','indexed_real_chm13','indexed_12_windows'):
    truth=json.loads((repo/'docs/coverage_oracles_2026-10-29_v1'/name/'truth.json').read_text())
    expected={(g,c,s,e if s=='+' else b) for o in truth for g,c,b,e,s in [o['key']]}
    lines=(root/(name+'.tsv')).read_text().splitlines(); i=next(i for i,line in enumerate(lines) if line.startswith('guide\tcontig\t'))
    rows=list(csv.DictReader(lines[i:],delimiter='\t')); got={(r['guide'],r['contig'],r['strand'],int(r['anchor'])) for r in rows}
    assert len(rows)==len(got) and got==expected, (name,len(expected-got),len(got-expected))
    print(name, 'expected/recovered', len(expected), 'duplicates',len(rows)-len(got))
lines=(root/'alt/out.tsv').read_text().splitlines(); i=next(i for i,line in enumerate(lines) if line.startswith('guide\tcontig\t'))
rows=list(csv.DictReader(lines[i:],delimiter='\t'))
for guide,contig,strand,start,end in [('g','terminal','+',2,20),('near','near','+',2,20),('g','reverse','-',3,21)]:
    matching=[r for r in rows if (r['guide'],r['contig'],r['strand'],int(r['start']),int(r['end']))==(guide,contig,strand,start,end)]
    assert matching and int(matching[0]['mismatches'])==0 and int(matching[0]['inserted_guide_bases'])==2
print('alternative alignment examples: 3/3')
for base in 'ACGT':
    for strand, suffix, anchor in [('+','plus',20),('-','minus',3)]:
        assert any(r['guide']=='g' and r['contig']=='pam_'+base+'_'+suffix and r['strand']==strand and int(r['anchor'])==anchor and r['pam']==base+'GG' for r in rows), (base,strand)
print('PAM wildcard bases on both strands: 8/8')
