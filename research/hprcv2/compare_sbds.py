#!/usr/bin/env python3
"""Read-only, bounded historical comparison on the completed HG00735 pilot.

Invalid old PAF geometry is labelled, never rewritten. Historical CIGAR witnesses
are diagnosed separately from standards-valid exports and from genome-wide recall.
"""
import argparse
import collections
import csv
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile

from collect import collect, only, site_records
from inventory import digest, library, SBDS


def replay(guide, target, cigar):
    if not re.fullmatch(r'(?:[1-9][0-9]*[=XID])+', cigar):
        raise ValueError('unsupported CIGAR')
    ops = [(int(n), op) for n, op in re.findall(r'(\d+)([=XID])', cigar)]
    q = r = matches = mismatches = gap_events = max_gap = 0
    aq, at = [], []
    for n, op in ops:
        if op in 'ID':
            gap_events += 1
            max_gap = max(max_gap, n)
        for _ in range(n):
            g = '-' if op == 'D' else guide[q]
            t = '-' if op == 'I' else target[r]
            if op in '=X':
                if (g == t) != (op == '='):
                    raise ValueError('CIGAR/FASTA assertion mismatch')
                matches += op == '='
                mismatches += op == 'X'
            aq.append(g)
            at.append(t)
            q += op != 'D'
            r += op != 'I'
    if q != len(guide) or r != len(target):
        raise ValueError('CIGAR consumption mismatch')
    return ops, aq, at, matches, mismatches, gap_events, max_gap


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    report = collect(a.root, 'sbds')
    if report['status'] != 'COMPLETE_SAFETY_PILOT':
        raise ValueError('both corrected SBDS pilot domains must be complete')
    plan = json.loads((a.root/'00_manifest/pilot-plan.json').read_text())
    base = a.root.parents[1]
    old = base/'results/SBDS_all_HPRCv2/HG00735_pat'
    old_guides = dict(library(old/'SBDS.fa'))
    if set(old_guides.values()) != {SBDS}:
        raise ValueError('historical guide identity mismatch')
    prior_pam = Path('/moosefs/salehi/SBDS_CRISPRapido/analysis/SBDS_validation_audit/per_assembly/HG00735_pat/PAM_audit_records.tsv')
    prior = {int(r['paf_line']): r for r in csv.DictReader(prior_pam.open(), delimiter='\t')}
    inventory = json.loads((a.root/'00_manifest/audit-v6/historical_sbds.json').read_text())['assemblies']['HG00735_pat']
    if digest(old/'SBDS.paf') != inventory['sha256']:
        raise ValueError('historical PAF changed since initial audit')
    sys.path.insert(0, str(Path(plan['source'])/'scripts/release_validation'))
    import audit_emitted as independent
    mm, pm = independent.matrices(plan['source'])
    lengths = {f[0]: int(f[1]) for f in (line.split('\t') for line in Path(plan['fai']).read_text().splitlines())}
    candidates = {}
    for task, completed in report['configurations'].items():
        marker = json.loads(Path(completed['completion_marker']).read_text())
        candidates[task] = site_records(only(marker['outputs'], 'eligible.both.paf'))
    a.out.mkdir()
    (a.out/'INCOMPLETE').write_text('Historical diagnostic comparison not committed.\n')
    counts, rows = collections.Counter(), []
    with tempfile.TemporaryDirectory(prefix='legacy-reference-', dir=str(a.out)) as temp:
        fasta = Path(temp)/'reference.fa'
        fasta.symlink_to(plan['reference'])
        shutil.copyfile(plan['fai'], str(fasta)+'.fai')
        if digest(plan['fai']) != plan['pins'][plan['fai']]:
            raise ValueError('FAI changed')
        if digest(plan['reference']) != plan['reference_sha256']:
            raise ValueError('reference changed')
        raw_rows = (old/'SBDS.paf').read_text().splitlines()
        for offset in range(0, len(raw_rows), 200):
            batch = []
            for n, raw in enumerate(raw_rows[offset:offset+200], offset+1):
                f = raw.split('\t')
                start, end = int(f[7]), int(f[8])
                target_region = '{}:{}-{}'.format(f[5], start+1, end)
                pam_region = ('{}:{}-{}'.format(f[5], end+1, end+3) if f[4] == '+' else '{}:{}-{}'.format(f[5], start-2, start))
                batch.append((n, raw, f, target_region, pam_region))
            fetched = independent.fetch_regions(plan['samtools'], fasta, [region for _, _, _, t, p in batch for region in (t, p)])
            for n, raw, f, target_region, pam_region in batch:
                target, pam = fetched[target_region], fetched[pam_region]
                if f[4] == '-':
                    target, pam = independent.rc(target), independent.rc(pam)
                tags = {s.split(':', 2)[0]: s.split(':', 2)[2] for s in f[12:]}
                cached = prior[n]
                if (cached['contig'], cached['ref_start0'], cached['ref_end0'], cached['strand'], cached['CIGAR']) != (f[5], f[7], f[8], f[4], tags['cg']):
                    raise ValueError('prior PAM audit does not match original PAF identity')
                if cached['target_guide_orientation'] != target or cached['PAM_3bp'] != pam:
                    raise ValueError('current FASTA differs from prior PAM audit')
                replay_error = ''
                try:
                    ops, aq, at, matched, mismatched, gaps, max_gap = replay(old_guides[f[0]], target, tags['cg'])
                except (ValueError, IndexError) as exc:
                    replay_error = str(exc)
                eligible = not replay_error and len(pam) == 3 and set(pam) <= set('ACGT') and pam[1:] == 'GG' and mismatched <= 2 and gaps <= 1 and max_gap <= 1 and matched/20 >= .75
                counts['old_cigar_replay_failures'] += bool(replay_error)
                try:
                    independent.parse_row(raw, old_guides, lengths)
                    geometry, geometry_reason = True, ''
                except independent.AuditError as exc:
                    geometry, geometry_reason = False, str(exc)
                key = ('SBDS', f[5], f[4], f[8] if f[4] == '+' else f[7])  # explicit exact-sequence mapping SBDS_therapeutic -> SBDS
                local = independent.independent_cfd(aq, at, ops, pam, mm, pm) if eligible else None
                row = dict(old_line=n, contig=f[5], start0=f[7], end0=f[8], strand=f[4],
                           target=target, pam=pam, old_cigar=tags['cg'], old_cfd=tags.get('cf', 'N/A'),
                           strict_m2_eligible_cigar_witness=eligible, standard_paf_geometry=geometry,
                           geometry_reason=geometry_reason, cigar_replay_error=replay_error,
                           independently_reconstructed_cfd=local if local is not None else 'N/A')
                counts['old_raw_rows'] += 1
                counts['old_strict_eligible_witnesses'] += eligible
                counts['old_geometry_failures'] += not geometry
                counts['old_canonical_NGG'] += pam[1:] == 'GG'
                counts['old_non_NGG'] += pam[1:] != 'GG'
                counts['eligible_old_rows_current_formula_unavailable'] += eligible and local is None
                counts['eligible_old_cfd_reconstruction_disagreements'] += eligible and local is not None and ('cf' not in tags or abs(float(tags['cf'])-local) > .000051)
                for task, records in candidates.items():
                    new = records.get(key)
                    row[task+'_same_anchor'] = new is not None
                    row[task+'_cigar'] = new['cigar'] if new else 'N/A'
                    row[task+'_cfd'] = new['cfd'] if new and new['cfd'] is not None else 'N/A'
                    counts[task+'_old_eligible_witnesses_missing'] += eligible and new is None
                    counts[task+'_raw_shared_anchors'] += new is not None
                    counts[task+'_shared_representative_differences'] += new is not None and (new['cigar'] != tags['cg'] or new['span'] != f[7:9])
                    counts[task+'_shared_stored_cfd_differences'] += new is not None and new['cfd'] is not None and 'cf' in tags and float(new['cfd']) != float(tags['cf'])
                rows.append(row)
        if digest(plan['reference']) != plan['reference_sha256'] or digest(plan['fai']) != plan['pins'][plan['fai']]:
            raise ValueError('reference/FAI changed during comparison')
    with (a.out/'records.tsv').open('x') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    summary = dict(status='PASS' if not any(v for k, v in counts.items() if k.endswith('_old_eligible_witnesses_missing')) else 'INVESTIGATE_MISSING',
                   counts=dict(counts), source_sha256={str(p): digest(p) for p in [old/'SBDS.paf', old/'SBDS.fa', prior_pam]},
                   reference_sha256=plan['reference_sha256'],
                   scope='One verified assembly. Legacy CIGAR witnesses separate from unmodified PAF geometry. Exact sequence-mapped guide and PAM-boundary keys, no crop/realignment or score-based representative selection. Not genome-wide recall.',
                   records_sha256=digest(a.out/'records.tsv'))
    with (a.out/'summary.json').open('x') as f:
        json.dump(summary, f, indent=2)
        f.write('\n')
    (a.out/'INCOMPLETE').unlink()
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
