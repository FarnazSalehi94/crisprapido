#!/usr/bin/env python3
"""Independent streaming PAF audit using indexed FASTA, not Rust target/PAM fields.

This validates emitted-row soundness. It cannot prove that a search found every
possible locus. Negative-strand CIGARs in this local release use guide orientation.
Outputs are exclusive; failed rows remain in errors.jsonl and no PASS is written.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys


class AuditError(ValueError):
    pass


def require(ok, message):
    if not ok:
        raise AuditError(message)


def rc(seq):
    return seq.translate(str.maketrans('ACGTRYMKSWBDHVN', 'TGCAYRKMSWVHDBN'))[::-1]


def load_guides(path):
    path = Path(path)
    if path.suffix == '.tsv':
        pairs = [(r['guide_id'], r['sequence']) for r in csv.DictReader(path.open(), delimiter='\t')]
    else:
        pairs, name, sequence = [], None, ''
        for line in path.read_text().splitlines():
            if line.startswith('>'):
                if name is not None:
                    pairs.append((name, sequence))
                name, sequence = line[1:].split()[0], ''
            else:
                require(name is not None, 'sequence before FASTA header')
                sequence += line.strip()
        if name is not None:
            pairs.append((name, sequence))
    require(bool(pairs), 'empty guide set')
    require(len({name for name, _ in pairs}) == len(pairs), 'duplicate guide ID')
    require(all(len(seq) == 20 and set(seq) <= set('ACGT') for _, seq in pairs), 'guides must be resolved 20nt DNA')
    return dict(pairs)


def matrices(directory):
    def read(name):
        pairs = [line.split() for line in (Path(directory)/name).read_text().splitlines() if line.strip()]
        require(all(len(p) == 2 for p in pairs), 'invalid matrix row')
        require(len({p[0] for p in pairs}) == len(pairs), 'duplicate matrix key')
        return {k: float(v) for k, v in pairs}
    return read('mismatch_scores.txt'), read('pam_scores.txt')


def parse_row(line, guides, lengths):
    f = line.rstrip('\r\n').split('\t')
    require(len(f) >= 12, 'PAF requires twelve fields')
    require(f[0] in guides, 'unknown guide ID')
    for i in [1, 2, 3, 6, 7, 8, 9, 10, 11]:
        require(bool(re.fullmatch(r'[0-9]+', f[i])), 'invalid numeric field '+str(i+1))
    require(f[1:4] == ['20', '0', '20'], 'invalid full-guide query coordinates')
    require(f[4] in ['+', '-'], 'invalid strand')
    require(f[5] in lengths and int(f[6]) == lengths[f[5]], 'contig ID/length mismatch')
    start, end = int(f[7]), int(f[8])
    require(0 <= start < end <= lengths[f[5]], 'invalid target interval')
    require(int(f[11]) == 255, 'unexpected mapping quality')
    require((end+3 <= lengths[f[5]]) if f[4] == '+' else (start >= 3), 'incomplete PAM interval')
    tags = {}
    for text in f[12:]:
        p = text.split(':', 2)
        require(len(p) == 3 and p[0] not in tags, 'duplicate or malformed auxiliary tag')
        tags[p[0]] = (p[1], p[2])
    require(tags.get('cg', ('', ''))[0] == 'Z', 'missing/wrong CIGAR tag')
    cigar = tags['cg'][1]
    require(bool(re.fullmatch(r'(?:[1-9][0-9]*[=XID])+', cigar)), 'unsupported or malformed CIGAR')
    ops = [(int(n), op) for n, op in re.findall(r'([0-9]+)([=XID])', cigar)]
    require(all(n <= 100 for n, _ in ops), 'implausibly large operation')
    require(sum(n for n, op in ops if op != 'D') == 20, 'CIGAR query consumption mismatch')
    require(sum(n for n, op in ops if op != 'I') == end-start, 'CIGAR reference consumption mismatch')
    target_region = '{}:{}-{}'.format(f[5], start+1, end)
    pam_region = ('{}:{}-{}'.format(f[5], end+1, end+3) if f[4] == '+' else
                  '{}:{}-{}'.format(f[5], start-2, start))
    return f, tags, ops, target_region, pam_region


def independent_cfd(aligned_guide, aligned_target, ops, pam, mm, pm):
    # Explicit existing CFD convention: a pure DNA bulge drops the same number
    # of PAM-distal columns as added DNA columns. No arbitrary right trimming.
    dna = sum(n for n, op in ops if op == 'D')
    if dna and not any(op == 'I' for _, op in ops):
        aligned_guide, aligned_target = aligned_guide[dna:], aligned_target[dna:]
    if len(aligned_guide) != 20 or pam[1:] not in pm:
        return None
    penalties = [pm[pam[1:]]]
    complement = dict(zip('ACGT-', 'TGCA-'))
    for position, (g, t) in enumerate(zip(aligned_guide, aligned_target), 1):
        if g == t or (position == 1 and '-' in (g, t)):
            continue
        key = 'r{}:d{},{}'.format(g.replace('T', 'U'), complement.get(t, t), position)
        if key not in mm:
            return None
        penalties.append(mm[key])
    product = 1.0
    for penalty in penalties:
        product *= penalty
    return product


def validate_alignment(parsed, guides, target, pam, mm, pm):
    f, tags, ops, _, _ = parsed
    if f[4] == '-':
        target, pam = rc(target), rc(pam)
    require(len(pam) == 3 and set(pam) <= set('ACGT') and pam[1:] == 'GG', 'noncanonical/ambiguous PAM')
    require(len(target) == int(f[8])-int(f[7]), 'FASTA target length mismatch')
    guide = guides[f[0]]
    q = r = 0
    aq, at = [], []
    matched = mismatched = gaps = max_gap = 0
    previous = None
    for length, op in ops:
        require(previous != op, 'unmerged adjacent CIGAR operations')
        previous = op
        if op in 'ID':
            gaps += 1
            max_gap = max(max_gap, length)
        for _ in range(length):
            if op in '=X':
                g, t = guide[q], target[r]
                require((g == t) if op == '=' else (g != t), 'CIGAR base assertion mismatch')
                q += 1
                r += 1
                matched += op == '='
                mismatched += op == 'X'
            elif op == 'I':
                g, t = guide[q], '-'
                q += 1
            else:
                g, t = '-', target[r]
                r += 1
            aq.append(g)
            at.append(t)
    require(q == 20 and r == len(target), 'sequence consumption mismatch')
    require(mismatched <= 4 and gaps <= 1 and max_gap <= 1 and matched >= 15, 'configured mismatch/bulge/match-fraction violation')
    require(int(f[9]) == matched and int(f[10]) == sum(n for n, _ in ops), 'PAF match/block count mismatch')
    expected_tags = dict(nm=mismatched, ng=gaps, bs=max_gap, **{'as': 3*mismatched+sum(5+n for n, op in ops if op in 'ID')})
    for name, value in expected_tags.items():
        require(tags.get(name) == ('i', str(value)), 'wrong alignment tag '+name)
    expected = independent_cfd(aq, at, ops, pam, mm, pm)
    cf = tags.get('cf')
    if cf is not None:
        require(cf[0] == 'f', 'wrong cf tag type')
        try:
            actual = float(cf[1])
        except ValueError:
            raise AuditError('invalid CFD number')
        require(math.isfinite(actual) and 0 <= actual <= 1, 'nonfinite/out-of-range CFD')
        require(expected is not None, 'unexpected score for unsupported reconstruction')
        require(abs(actual-expected) <= .000051, 'CFD score differs from independent reconstruction')
    else:
        actual = None
        require(expected is None, 'missing score for scoreable row')
    return dict(guide=f[0], strand=f[4], contig=f[5], start0=int(f[7]), end0=int(f[8]),
                pam=pam, target=target, aligned_guide=''.join(aq), aligned_target=''.join(at),
                mismatches=mismatched, bulges=gaps, cfd=actual, independent_cfd=expected,
                status='SCORED' if actual is not None else 'UNSUPPORTED_CFD')


def fetch_regions(samtools, fasta, regions):
    unique = sorted(set(regions))
    if not unique:
        return {}
    result = subprocess.run([samtools, 'faidx', str(fasta)]+unique,
                            capture_output=True, text=True)
    require(result.returncode == 0 and not result.stderr, 'faidx failed: '+result.stderr)
    sequences, name = {}, None
    for line in result.stdout.splitlines():
        if line.startswith('>'):
            name = line[1:]
            sequences[name] = ''
        else:
            require(name is not None, 'malformed faidx response')
            sequences[name] += line.strip().upper()
    require(set(sequences) == set(unique), 'missing faidx interval')
    return sequences


def audit(paf, fasta, guide_path, matrix_dir, samtools, output, batch_size=500):
    require(1 <= batch_size <= 1000, 'batch size must be 1..1000')
    output = Path(output)
    output.mkdir(exist_ok=False)
    guides, (mm, pm) = load_guides(guide_path), matrices(matrix_dir)
    lengths = {x[0]: int(x[1]) for x in (s.split('\t') for s in Path(str(fasta)+'.fai').read_text().splitlines())}
    counts = collections.Counter()
    per_guide = {name: collections.Counter() for name in guides}
    summary = dict(status='RUNNING', rows=0, errors=0)
    fields = 'line guide strand contig start0 end0 pam target aligned_guide aligned_target mismatches bulges cfd independent_cfd status'.split()
    with (output/'records.tsv').open('w') as records, (output/'errors.jsonl').open('w') as errors:
        writer = csv.DictWriter(records, delimiter='\t', fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        def error(number, line, reason):
            summary['errors'] += 1
            errors.write(json.dumps(dict(line=number, error=str(reason), paf=line.rstrip()))+'\n')
        def process(batch):
            parsed, regions = [], []
            for number, line in batch:
                try:
                    p = parse_row(line, guides, lengths)
                    parsed.append((number, line, p))
                    regions.extend(p[3:5])
                except (ValueError, IndexError) as e:
                    error(number, line, e)
            try:
                sequences = fetch_regions(samtools, fasta, regions)
            except AuditError as e:
                for number, line, _ in parsed:
                    error(number, line, e)
                return
            for number, line, p in parsed:
                try:
                    result = validate_alignment(p, guides, sequences[p[3]], sequences[p[4]], mm, pm)
                    writer.writerow(dict(line=number, **result))
                    counts['verified'] += 1
                    counts['strand_'+result['strand']] += 1
                    counts[result['status']] += 1
                    counts['bulged' if result['bulges'] else 'ungapped'] += 1
                    counts['CFD_ge_020'] += result['cfd'] is not None and result['cfd'] >= .20
                    counts['CFD_lt_020'] += result['cfd'] is not None and result['cfd'] < .20
                    per_guide[result['guide']]['verified'] += 1
                    per_guide[result['guide']]['strand_'+result['strand']] += 1
                except (ValueError, IndexError) as e:
                    error(number, line, e)
        try:
            with open(paf) as f:
                batch = []
                for number, line in enumerate(f, 1):
                    summary['rows'] += 1
                    batch.append((number, line))
                    if len(batch) == batch_size:
                        process(batch)
                        batch = []
                if batch:
                    process(batch)
        except (OSError, UnicodeError) as e:
            error(0, '', e)
    require(counts['verified']+summary['errors'] == summary['rows'], 'internal audit accounting error')
    def digest(path):
        h = hashlib.sha256()
        with Path(path).open('rb') as f:
            for block in iter(lambda: f.read(1024*1024), b''):
                h.update(block)
        return h.hexdigest()
    summary['provenance'] = dict(
        paf=str(Path(paf).resolve()), fasta=str(Path(fasta).resolve()),
        guides=str(Path(guide_path).resolve()), matrices=str(Path(matrix_dir).resolve()),
        samtools=str(samtools), batch_size=batch_size,
        auditor_sha256=digest(__file__), paf_sha256=digest(paf),
        guide_sha256=digest(guide_path), fai_sha256=digest(str(fasta)+'.fai'),
        matrix_sha256={name: digest(Path(matrix_dir)/name) for name in
                       ['mismatch_scores.txt', 'pam_scores.txt']})
    summary.update(counts=dict(counts), per_guide={k: dict(v) for k, v in per_guide.items()},
                   status='PASS' if summary['errors'] == 0 else 'FAIL',
                   limitation='Emitted-row soundness only; not search sensitivity/completeness.')
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    if summary['status'] == 'PASS':
        (output/'PASS').write_text('Every emitted row independently verified.\n')
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['paf', 'fasta', 'guides', 'matrices', 'samtools', 'output']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--batch-size', type=int, default=500)
    a = p.parse_args()
    result = audit(a.paf, a.fasta, a.guides, a.matrices, a.samtools, a.output, a.batch_size)
    print(json.dumps(result, indent=2))
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
