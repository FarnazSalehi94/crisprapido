#!/usr/bin/env python3
"""Record full evidence and compare loci anchored by guide, strand and exact PAM."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
from coverage_oracle import enumerate_sites, rc


def fasta(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            result[name] = ''
        else:
            result[name] += line.strip().upper()
    return result


def category(row):
    return (('exact' if row['mismatches'] == 0 else 'mismatch') if row['bulge'] == 'none'
            else ('I' if row['bulge'] == 'RNA' else 'D') if row['mismatches'] == 0
            else 'mismatch_bulge')


def locus(key, refs):
    guide, contig, start, end, strand = key
    bounds = (end, end + 3) if strand == '+' else (start - 3, start)
    seq = refs[contig][bounds[0]:bounds[1]] if bounds[0] >= 0 else ''
    pam = seq if strand == '+' else rc(seq)
    return (guide, contig, strand, bounds[0], bounds[1], pam)


def ops_for_alignment(a, b):
    out = []
    for x, y in zip(a, b):
        op = 'D' if x == '-' else 'I' if y == '-' else '=' if x == y else 'X'
        if out and out[-1][1] == op:
            out[-1][0] += 1
        else:
            out.append([1, op])
    return ''.join(str(n) + op for n, op in out)


def replay(guide, target, cigar, reverse_sam=False):
    parts = [(int(n), op) for n, op in re.findall(r'(\d+)([MID=X])', cigar)]
    if ''.join(str(n)+op for n, op in parts) != cigar:
        return {'error': 'unsupported CIGAR: ' + cigar}
    if reverse_sam:
        parts.reverse()
    qi = ti = 0
    a, b = [], []
    for n, op in parts:
        for _ in range(n):
            if op in 'M=X':
                if qi >= len(guide) or ti >= len(target):
                    return {'error': 'CIGAR consumes outside guide or target'}
                a.append(guide[qi]); b.append(target[ti]); qi += 1; ti += 1
            elif op == 'I':
                if qi >= len(guide): return {'error': 'CIGAR consumes outside guide'}
                a.append(guide[qi]); b.append('-'); qi += 1
            else:
                if ti >= len(target): return {'error': 'CIGAR consumes outside target'}
                a.append('-'); b.append(target[ti]); ti += 1
    if qi != len(guide) or ti != len(target):
        return {'error': 'CIGAR does not consume full guide and target'}
    a, b = ''.join(a), ''.join(b)
    return dict(aligned_guide=a, aligned_target=b, normalized_cigar=ops_for_alignment(a, b),
                mismatches=sum(x != y and x != '-' and y != '-' for x, y in zip(a, b)),
                gap_positions=[i for i, (x, y) in enumerate(zip(a, b)) if x == '-' or y == '-'])


def identity(key, refs):
    guide, contig, start, end, strand = key
    target = refs[contig][start:end]
    pam_key = locus(key, refs)
    return dict(guide=guide, contig=contig, start=start, end=end, strand=strand,
                reference_span=end-start, reference_sequence=target,
                oriented_target=target if strand == '+' else rc(target),
                pam_start=pam_key[3], pam_end=pam_key[4], pam=pam_key[5])


def overlaps(a, b):
    return a[0] == b[0] and a[1] == b[1] and a[4] == b[4] and a[2] < b[3] and b[2] < a[3]


def analyze(directory):
    d = Path(directory)
    refs, guides = fasta(d/'ref.fa'), fasta(d/'guides.fa')
    oracle = defaultdict(list)
    for guide, sequence in guides.items():
        for contig, ref in refs.items():
            for row in enumerate_sites(ref, sequence):
                key = (guide, contig, row['start'], row['end'], row['strand'])
                row['cigar'] = ops_for_alignment(row['aligned_guide'], row['aligned_target'])
                oracle[(key, category(row))].append(row)
    final = []
    for line in (d/'out.paf').read_text().splitlines():
        f = line.split('\t'); key = (f[0], f[5], int(f[7]), int(f[8]), f[4])
        tags = {x[:2]: x[5:] for x in f[12:]}
        details = identity(key, refs)
        details.update(cigar=tags['cg'], alignment=replay(guides[key[0]], details['oriented_target'], tags['cg']))
        if 'error' in details['alignment']: raise ValueError(details)
        final.append((key, details))
    raw = []
    for line in (d/'edit.sam').read_text().splitlines():
        if line.startswith('@'): continue
        f = line.split('\t'); flag = int(f[1])
        if flag & 4: continue
        span = sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]) if op in 'MDN=X')
        key = (f[0], f[2], int(f[3])-1, int(f[3])-1+span, '-' if flag & 16 else '+')
        details = identity(key, refs)
        details.update(cigar_sam=f[5], tags=f[11:],
                       alignment=replay(guides[key[0]], details['oriented_target'], f[5], key[4]=='-'))
        raw.append((key, details))
    existing = json.loads((d/'discrepancies.json').read_text())
    rows = []
    for miss in existing:
        key = tuple(miss['key']); cls = miss['category']
        candidates = [v for k, v in raw if overlaps(key, k)]
        reports = [v for k, v in final if overlaps(key, k)]
        same = [(k, v) for k, v in final if locus(k, refs) == locus(key, refs)]
        # B is a different representation at the same PAM-anchored site,
        # not a claim of identical edit structure or biological outcome.
        classification = 'B' if same else 'C' if reports else 'E'
        record = identity(key, refs)
        record.update(category=cls, oracle_alignments=oracle[(key, cls)],
                      overlapping_raw=candidates, overlapping_final=reports,
                      same_pam_reported=[v for _, v in same], classification=classification)
        if not record['oracle_alignments']: raise ValueError('missing oracle representation')
        rows.append(record)
    truth_loci = {locus(k, refs) for k, _ in oracle}
    final_loci = {locus(k, refs) for k, _ in final}
    metrics = {'oracle_loci':len(truth_loci), 'reported_loci':len(final_loci),
               'recovered_loci':len(truth_loci & final_loci),
               'recall':len(truth_loci & final_loci)/len(truth_loci) if truth_loci else None,
               'missing_loci':sorted(truth_loci-final_loci),
               'unexplained_final_loci':sorted(final_loci-truth_loci),
               'class_counts':dict(Counter(x['classification'] for x in rows)),
               'raw_replay_errors':sum('error' in x['alignment'] for _, x in raw)}
    return rows, metrics


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory'); p.add_argument('--output', required=True)
    a = p.parse_args()
    rows, metrics = analyze(a.directory)
    Path(a.output).write_text(json.dumps({'metrics':metrics, 'discrepancies':rows}, indent=2)+'\n')
    print(json.dumps(metrics, sort_keys=True))

if __name__ == '__main__': main()
