#!/usr/bin/env python3
"""Small, indexed-FASTA regression comparison. Never scans an assembly or overwrites output.

Requires an audited record table, guide table, exact indexed source FASTA and
old/new executables. The old executable may need an explicit loader prefix.
See docs/PAM_CFD_VALIDATION.md for the recorded validation run.
"""
import argparse
import collections
import csv
import hashlib
import json
import pathlib
import random
import re
import subprocess


def rc(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def fetch(samtools, fasta, regions):
    result = subprocess.run([samtools, 'faidx', str(fasta)] + list(dict.fromkeys(regions)),
                            text=True, capture_output=True, check=True)
    if result.stderr:
        raise RuntimeError(result.stderr)
    seqs = {}
    for line in result.stdout.splitlines():
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = ''
        else:
            seqs[name] += line.strip().upper()
    return seqs


def score_row(f, guide, target, pam, mm_scores, pam_scores):
    cigar = next(t[5:] for t in f[12:] if t.startswith('cg:Z:'))
    q = r = 0
    aligned_q, aligned_r = [], []
    dna = 0
    has_i = False
    for n, op in re.findall(r'(\d+)([=XIDM])', cigar):
        for _ in range(int(n)):
            if op in '=XM':
                a, b = guide[q], target[r]
                if op == '=':
                    assert a == b, (f, q, r, a, b)
                if op == 'X':
                    assert a != b, (f, q, r, a, b)
                q += 1
                r += 1
            elif op == 'I':
                a, b = guide[q], '-'
                q += 1
                has_i = True
            else:
                a, b = '-', target[r]
                r += 1
                dna += 1
            aligned_q.append(a)
            aligned_r.append(b)
    assert q == len(guide) and r == len(target)
    assert f[2:4] == ['0', str(len(guide))]
    if dna and not has_i:
        aligned_q = aligned_q[dna:]
        aligned_r = aligned_r[dna:]
    if len(aligned_q) != 20 or len(pam) != 3 or set(pam) - set('ACGT'):
        return None
    score = pam_scores[pam[1:]]
    for i, (a, b) in enumerate(zip(aligned_q, aligned_r)):
        if a == b or (i == 0 and '-' in (a, b)):
            continue
        key = 'r{}:d{},{}'.format(a.replace('T', 'U'), rc(b), i + 1)
        if key not in mm_scores:
            return None
        score *= mm_scores[key]
    return score


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ['records', 'guides', 'assembly', 'fasta', 'old-bin', 'new-bin', 'samtools', 'matrices', 'output']:
        p.add_argument('--' + arg, required=True)
    p.add_argument('--old-prefix', action='append', default=[], help='Repeat for each loader argument')
    p.add_argument('--limit', type=int, default=120)
    a = p.parse_args()
    if not 1 <= a.limit <= 1000:
        p.error('limit must be between 1 and 1000; this is a small validation only')
    out = pathlib.Path(a.output)
    out.mkdir(exist_ok=False)
    fa = pathlib.Path(a.fasta)
    assert pathlib.Path(str(fa) + '.fai').exists()
    assert fa.suffix != '.gz' or pathlib.Path(str(fa) + '.gzi').exists()
    rows = [r for r in csv.DictReader(open(a.records), delimiter='\t') if r['assembly'] == a.assembly]
    rng = random.Random(20261028)
    selected = []
    # Preserve every alignment category, then fill uniformly without count weighting.
    groups = collections.defaultdict(list)
    for row in rows:
        c = row['cigar']
        group = (row['strand'], 'leading_D' if c.startswith('1D') else
                 'trailing_D' if c.endswith('1D') else 'DNA' if 'D' in c else
                 'RNA' if 'I' in c else 'mismatch' if 'X' in c else 'perfect')
        groups[group].append(row)
    for group in sorted(groups):
        selected.extend(rng.sample(groups[group], min(5, len(groups[group]))))
    rest = [r for r in rows if r not in selected]
    selected = (selected + rng.sample(rest, min(len(rest), max(0, a.limit-len(selected)))))[:a.limit]
    guides = {r['guide_id']: r['sequence'] for r in csv.DictReader(open(a.guides), delimiter='\t')}
    lengths = {r.split('\t')[0]: int(r.split('\t')[1]) for r in pathlib.Path(str(fa)+'.fai').read_text().splitlines()}
    regions = []
    coordinates = []
    for row in selected:
        # Historical reporting advanced start for a leading reference deletion.
        leading = re.match(r'^(\d+)D', row['cigar'])
        shift = int(leading.group(1)) if leading else 0
        s, e = int(row['tstart0'])-shift, int(row['tend0'])-shift
        left, right = max(0, s-6), min(lengths[row['contig']], e+6)
        regions.append('{}:{}-{}'.format(row['contig'], left+1, right))
        coordinates.append((s-left, e-left))
    seqs = fetch(a.samtools, fa, regions)
    mini, sam = [], []
    for i, (row, region, (s, e)) in enumerate(zip(selected, regions, coordinates)):
        name = 'row_{:04d}'.format(i)
        mini += ['>'+name, seqs[region]]
        sam.append('\t'.join(map(str, [row['guide_id'], 0 if row['strand']=='+' else 16, name, s+1, 255,
                                         row['cigar'], '*', 0, 0, '*', '*', 'NM:i:0', 'AS:i:0'])))
        row['mini_reference'] = name
        row['source_region'] = region
    mini_path = out/'subset.fa'
    mini_path.write_text('\n'.join(mini)+'\n')
    (out/'subset.sam').write_text('\n'.join(sam)+'\n')
    (out/'guides.fa').write_text(''.join('>{}\n{}\n'.format(k, guides[k]) for k in sorted({r['guide_id'] for r in selected})))
    with (out/'selected_records.tsv').open('w') as h:
        w = csv.DictWriter(h, fieldnames=list(selected[0]), delimiter='\t'); w.writeheader(); w.writerows(selected)
    subprocess.run([a.samtools, 'faidx', str(mini_path)], check=True)
    mm = dict((k, float(v)) for k, v in (l.split() for l in (pathlib.Path(a.matrices)/'mismatch_scores.txt').read_text().splitlines()))
    pm = dict((k, float(v)) for k, v in (l.split() for l in (pathlib.Path(a.matrices)/'pam_scores.txt').read_text().splitlines()))
    summary = {'assembly': a.assembly, 'input_records': len(selected), 'seed': 20261028, 'runs': {}}
    outputs = {}
    verification_rows = []
    for mode, binary, prefix, extra in [('old', a.old_bin, a.old_prefix, []), ('corrected_strict', a.new_bin, [], []),
                                       ('corrected_permissive', a.new_bin, [], ['--allow-noncanonical-pam'])]:
        cmd = prefix + [binary, '-r', str(mini_path), '-G', str(out/'guides.fa'), '--columba-sam', str(out/'subset.sam'),
                        '-p', 'GG', '-m', '4', '-b', '1', '-z', '1', '-t', '1',
                        '--mismatch-scores', str(pathlib.Path(a.matrices)/'mismatch_scores.txt'),
                        '--pam-scores', str(pathlib.Path(a.matrices)/'pam_scores.txt')] + extra
        result = subprocess.run(cmd, text=True, capture_output=True, check=True)
        (out/(mode+'.paf')).write_text(result.stdout)
        (out/(mode+'.stderr')).write_text(result.stderr)
        (out/(mode+'.command.json')).write_text(json.dumps(cmd, indent=2))
        records = [l.split('\t') for l in result.stdout.splitlines()]
        canonical = scored = retained = verified = 0
        needed = []
        for f in records:
            s, e = int(f[7]), int(f[8])
            name = f[5]
            needed += ['{}:{}-{}'.format(name,s+1,e)]
            needed += ['{}:{}-{}'.format(name,e+1,e+3) if f[4]=='+' else '{}:{}-{}'.format(name,s-2,s)]
        fetched = fetch(a.samtools, mini_path, needed) if needed else {}
        for i, f in enumerate(records):
            target, pam = fetched[needed[2*i]], fetched[needed[2*i+1]]
            if f[4]=='-':
                target, pam = rc(target), rc(pam)
            canonical += len(pam)==3 and pam[1:]=='GG'
            cf = next((float(t[5:]) for t in f[12:] if t.startswith('cf:f:')), None)
            scored += cf is not None
            retained += cf is not None and cf >= .20
            if mode != 'old':
                expected = score_row(f, guides[f[0]], target, pam, mm, pm)
                assert (expected is None) == (cf is None), (f, expected)
                if expected is not None:
                    assert abs(cf-expected) <= .000051, (f,expected)
                if mode == 'corrected_strict':
                    assert pam[1:] == 'GG'
                verification_rows.append({'mode': mode, 'guide': f[0], 'reference': f[5], 'strand': f[4],
                                          'start0': f[7], 'end0': f[8], 'target': target, 'pam': pam,
                                          'reported_CFD': cf, 'independent_CFD': expected})
                verified += 1
        outputs[mode] = records
        summary['runs'][mode] = {'records': len(records), 'canonical_NGG': canonical, 'non_NGG': len(records)-canonical,
                                'scored': scored, 'CFD_ge_020': retained, 'independently_verified': verified}
    def identity(f): return tuple(f[:1]+f[4:6]+f[7:9])
    old = {identity(f): f for f in outputs['old']}
    for mode in ['corrected_strict','corrected_permissive']:
        shared = changed = 0
        for f in outputs[mode]:
            if identity(f) in old:
                shared += 1
                cf = lambda row: next((t for t in row[12:] if t.startswith('cf:f:')), None)
                changed += cf(f) != cf(old[identity(f)])
        summary['runs'][mode]['shared_coordinate_records_with_old'] = shared
        summary['runs'][mode]['shared_records_changed_CFD'] = changed
    for mode in ['corrected_strict', 'corrected_permissive']:
        by_ref = lambda records: {f[5]: f for f in records}
        old_refs, new_refs = by_ref(outputs['old']), by_ref(outputs[mode])
        changes = collections.Counter()
        get_cf = lambda f: next((float(t[5:]) for t in f[12:] if t.startswith('cf:f:')), None)
        for name in old_refs.keys() & new_refs.keys():
            old, new = old_refs[name], new_refs[name]
            before, after = get_cf(old), get_cf(new)
            changes['interval_changed'] += old[7:9] != new[7:9]
            if before is not None and after is not None:
                changes['score_increased'] += after > before
                changes['score_decreased'] += after < before
                changes['score_unchanged'] += after == before
        changes['removed_candidate_ids'] = len(old_refs.keys() - new_refs.keys())
        summary['runs'][mode]['by_candidate_id'] = dict(changes)
    if verification_rows:
        with (out/'independent_verification.tsv').open('w') as h:
            w = csv.DictWriter(h, fieldnames=list(verification_rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(verification_rows)
    summary['binary_sha256'] = {mode: hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
                               for mode,path in [('old',a.old_bin),('corrected',a.new_bin)]}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
