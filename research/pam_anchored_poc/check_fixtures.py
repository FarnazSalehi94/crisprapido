#!/usr/bin/env python3
"""Check all saved difficult panel anchor keys and existing fallback examples."""
import csv, json, re, sys
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
def read_fasta(path):
    out={};name=None
    for line in path.read_text().splitlines():
        if line.startswith('>'):name=line[1:].split()[0];out[name]=''
        else:out[name]+=line.upper()
    return out
refs=read_fasta(root/'alt/ref.fa');guides=read_fasta(root/'alt/guides.fa')
for guide,contig,strand,start,end in [('g','terminal','+',2,20),('near','near','+',2,20),('g','reverse','-',3,21)]:
    anchor=end if strand=='+' else start
    matching=[r for r in rows if (r['guide'],r['contig'],r['strand'],int(r['anchor']))==(guide,contig,strand,anchor)]
    assert matching and int(matching[0]['mismatches'])==0 and int(matching[0]['match_count'])>=15 and int(matching[0]['bulge_events'])<=1
    # Replay the chosen alternative; it need not retain the original 2I span.
    r=matching[0];seq=refs[contig][int(r['start']):int(r['end'])]
    if strand=='-':seq=seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
    q=guides[guide];qi=ti=matches=events=gaplen=0;prev=''
    ops=re.findall(r'(\d+)([=XID])',r['cigar']);assert ''.join(n+op for n,op in ops)==r['cigar']
    for raw,op in ops:
        n=int(raw)
        if op in 'ID':
            if not prev or prev not in 'ID':events+=1;gaplen=0
            gaplen+=n;assert gaplen<=2
        else:gaplen=0
        prev=op
        for _ in range(n):
            if op=='I':qi+=1
            elif op=='D':ti+=1
            else:assert op=='=' and q[qi]==seq[ti];qi+=1;ti+=1;matches+=1
    assert qi==len(q) and ti==len(seq) and events<=1 and matches>=15
print('alternative alignment examples: 3/3')
for base in 'ACGT':
    for strand, suffix, anchor in [('+','plus',20),('-','minus',3)]:
        assert any(r['guide']=='g' and r['contig']=='pam_'+base+'_'+suffix and r['strand']==strand and int(r['anchor'])==anchor and r['pam']==base+'GG' for r in rows), (base,strand)
print('PAM wildcard bases on both strands: 8/8')
