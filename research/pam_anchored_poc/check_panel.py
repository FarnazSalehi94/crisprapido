#!/usr/bin/env python3
"""Independently audit prototype rows and compare anchor sets with a retained oracle."""
import argparse, csv, json, re
from collections import Counter
from pathlib import Path

def fasta(path):
    out = {}; name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]; out[name] = ''
        else: out[name] += line.strip().upper()
    return out

def rc(seq): return seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('panel'); ap.add_argument('output'); args=ap.parse_args()
    p=Path(args.panel); refs=fasta(p/'ref.fa'); guides=fasta(p/'guides.fa')
    rows=Path(args.output).read_text().splitlines(); header=next(i for i,s in enumerate(rows) if s.startswith('guide\tcontig\t'))
    hits=list(csv.DictReader(rows[header:],delimiter='\t')); keys=[]
    for h in hits:
        g=guides[h['guide']]; seq=refs[h['contig']]; start,end=map(int,(h['start'],h['end'])); strand=h['strand']; anchor=int(h['anchor'])
        assert 0<=start<end<=len(seq)
        if strand=='+': target=seq[start:end]; pam=seq[end:end+3]; assert anchor==end
        else: target=rc(seq[start:end]); pam=rc(seq[start-3:start]); assert start>=3 and anchor==start
        assert len(pam)==3 and pam[0] in 'ACGT' and pam[1:]=='GG' and pam==h['pam']
        ops=re.findall(r'(\d+)([=XID])',h['cigar']); assert ''.join(n+op for n,op in ops)==h['cigar']
        gi=ti=mm=matches=ins=dele=groups=0
        for n,op in ops:
            n=int(n)
            if op in 'ID': groups+=1; assert n<=2
            for _ in range(n):
                if op=='I': gi+=1; ins+=1
                elif op=='D': ti+=1; dele+=1
                else:
                    assert gi<len(g) and ti<len(target) and ((g[gi]==target[ti])==(op=='=')); gi+=1;ti+=1; mm+=op=='X'; matches+=op=='='
        assert gi==len(g) and ti==len(target) and mm<=5 and groups<=1 and matches>=0.75*len(g)
        assert (mm,ins,dele)==tuple(int(h[k]) for k in ('mismatches','inserted_guide_bases','deleted_guide_bases'))
        keys.append((h['guide'],h['contig'],strand,anchor-(3 if strand=='-' else 0)))
    expected={(x['guide'],x['contig'],x['strand'],x['anchor']) for x in map(json.loads,(p/'oracle.jsonl').open())}
    actual=set(keys); assert len(keys)==len(actual), 'duplicate anchors'; assert actual==expected, (len(expected-actual),len(actual-expected))
    print(json.dumps({'expected':len(expected),'reported':len(keys),'recovered':len(expected&actual),'unexpected':len(actual-expected),'duplicates':len(keys)-len(actual),'strand_counts':dict(Counter(k[2] for k in keys))},sort_keys=True))
if __name__=='__main__': main()
