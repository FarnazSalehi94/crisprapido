#!/usr/bin/env python3
"""Independent small-window strict NGG coverage enumeration; coordinates zero-based half-open."""
import argparse
import json
import subprocess

DNA = set('ACGT')
TRANS = str.maketrans('ACGTN', 'TGCAN')

def rc(seq):
    return seq.upper().translate(TRANS)[::-1]

def enumerate_sites(sequence, guide, max_mismatch=4, bulges=True, min_fraction=.75):
    sequence, guide = sequence.upper(), guide.upper()
    if not guide or set(guide) - DNA or not 0 <= max_mismatch <= 4:
        raise ValueError('guide must be ACGT and mismatch limit in 0..4')
    n, L = len(sequence), len(guide)
    for strand in ('+', '-'):
        oriented = sequence if strand == '+' else rc(sequence)
        for start in range(n):
            for length in ((L-1, L, L+1) if bulges else (L,)):
                if length <= 0 or start + length + 3 > n:
                    continue
                target = oriented[start:start+length]
                pam = oriented[start+length:start+length+3]
                if set(target+pam)-DNA or pam[1:] != 'GG':
                    continue
                variants = []
                if length == L:
                    variants.append((guide, target, 'none', None))
                elif length == L-1:
                    variants.extend((guide, target[:i]+'-'+target[i:], 'RNA', i) for i in range(L))
                else:
                    variants.extend((guide[:i]+'-'+guide[i:], target, 'DNA', i) for i in range(L+1))
                for a,b,kind,pos in variants:
                    mismatches = sum(x != y for x,y in zip(a,b) if x != '-' and y != '-')
                    matches = sum(x == y for x,y in zip(a,b) if x != '-')
                    if mismatches <= max_mismatch and matches / L >= min_fraction:
                        s,e = (start,start+length) if strand == '+' else (n-start-length,n-start)
                        yield dict(start=s,end=e,strand=strand,pam=pam,mismatches=mismatches,bulge=kind,bulge_position=pos,aligned_guide=a,aligned_target=b)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',required=True); p.add_argument('--region',required=True,help='samtools 1-based inclusive region'); p.add_argument('--guide',required=True)
    p.add_argument('--samtools',default='samtools'); p.add_argument('--mode',choices=['exact','bounded'],default='bounded'); p.add_argument('--offset',type=int,required=True,help='zero-based window start; must match region')
    a=p.parse_args()
    chrom,span=a.region.rsplit(':',1); first,last=map(int,span.split('-'))
    if first-1 != a.offset or last < first: p.error('offset must match region start')
    out=subprocess.run([a.samtools,'faidx',a.fasta,a.region],check=True,capture_output=True,text=True).stdout
    seq=''.join(out.splitlines()[1:])
    if len(seq)!=last-first+1: raise ValueError('incomplete FASTA fetch')
    for row in enumerate_sites(seq,a.guide,0 if a.mode=='exact' else 4,a.mode!='exact'):
        row['start']+=a.offset; row['end']+=a.offset; row['contig']=chrom
        print(json.dumps(row,sort_keys=True))
if __name__=='__main__': main()
