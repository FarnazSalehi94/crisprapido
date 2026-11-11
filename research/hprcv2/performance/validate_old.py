#!/usr/bin/env python3
"""Label preserved historical PAF rows; count strict eligible witnesses and containment.

Rows are never rewritten or repaired. Geometry, replay, PAM and eligibility are
separate labels. This validator is a diagnostic for a preserved legacy binary whose
output is known to include non-NGG PAMs and legacy bulge query coordinates; it is
never a recall standard for the corrected implementation.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
for candidate in (HERE, HERE.parent):
    sys.path.insert(0, str(candidate))
from compare_sbds import replay  # noqa: E402
from inventory import digest, library  # noqa: E402

FIELDS = ['line', 'guide', 'contig', 'strand', 'start0', 'end0', 'qstart', 'qend', 'pam',
          'cigar', 'geometry', 'replay_error', 'mismatches', 'bulge_events', 'max_bulge_size',
          'strict_ngg', 'strict_eligible', 'stored_cfd', 'reconstructed_cfd', 'cfd_match', 'contained_in_new']


def label_row(number, fields, guides, target, pam, ops_cfd):
    """Pure labeling of one already-fetched historical row; no I/O, no repair."""
    result = dict.fromkeys(FIELDS[1:], '')
    if len(fields) < 12:
        result.update(line=number, geometry='paf_has_fewer_than_12_fields')
        return result, True
    result.update(line=number, guide=fields[0], contig=fields[5], strand=fields[4],
                  start0=fields[7], end0=fields[8], qstart=fields[2], qend=fields[3], pam=pam)
    tags = {}
    for item in fields[12:]:
        parts = item.split(':', 2)
        if len(parts) != 3 or parts[0] in tags:
            result['geometry'] = 'invalid_auxiliary_tag'
            return result, True
        tags[parts[0]] = parts[2]
    result['cigar'] = tags.get('cg', '')
    result['stored_cfd'] = tags.get('cf', 'N/A')
    if fields[1:4] == ['20', '0', '20']:
        result['geometry'] = 'standard'
    else:
        result['geometry'] = 'legacy_query_coordinates'
    try:
        int(fields[7]), int(fields[8])
    except ValueError:
        result['geometry'] += '+invalid_target_coordinates'
        return result, True
    result['strict_ngg'] = len(pam) == 3 and set(pam) <= set('ACGT') and pam[1:] == 'GG'
    if fields[0] not in guides or fields[4] not in '+-':
        result['replay_error'] = 'unknown_guide_or_strand'
        return result, True
    try:
        ops, aligned_guide, aligned_target, matched, mismatches, gaps, max_gap = replay(guides[fields[0]], target, result['cigar'])
    except (ValueError, IndexError) as exc:
        # A replay failure is an expected legacy diagnostic, never a validation abort.
        result['replay_error'] = str(exc)
        return result, True
    result.update(mismatches=mismatches, bulge_events=gaps, max_bulge_size=max_gap)
    eligible = (result['strict_ngg'] and mismatches <= 2 and gaps <= 1 and max_gap <= 1 and matched/20 >= .75)
    result['strict_eligible'] = eligible
    if eligible and result['stored_cfd'] != 'N/A' and ops_cfd is not None:
        expected = ops_cfd(aligned_guide, aligned_target, ops, pam)
        if expected is not None:
            result['reconstructed_cfd'] = expected
            result['cfd_match'] = abs(float(result['stored_cfd'])-expected) <= .000051
    return result, True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['paf', 'guides', 'reference', 'matrices', 'samtools', 'release-scripts', 'out']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--new-paf', required=True)
    a = p.parse_args()
    sys.path.insert(0, str(Path(a.release_scripts).resolve()))
    import audit_emitted as independent
    mm, pm = independent.matrices(a.matrices)
    guides = dict(library(a.guides))
    lengths = {f[0]: int(f[1]) for f in (line.split('\t') for line in Path(str(a.reference)+'.fai').read_text().splitlines())}
    new_keys = set()
    with open(a.new_paf) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            new_keys.add((fields[0], fields[5], fields[4], fields[8] if fields[4] == '+' else fields[7]))
    out = Path(a.out)
    out.mkdir()
    (out/'INCOMPLETE').write_text('Historical validation not committed.\n')

    def ops_cfd(aligned_guide, aligned_target, ops, pam):
        return independent.independent_cfd(aligned_guide, aligned_target, ops, pam, mm, pm)

    counts = dict(rows=0, accounted=0, strict_eligible=0, strict_eligible_contained_in_new=0,
                  strict_eligible_missing_in_new=0, rows_contained_in_new=0, cfd_compared=0,
                  cfd_disagreements=0)
    with open(a.paf) as source, (out/'records.tsv').open('x') as records:
        writer = csv.DictWriter(records, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        batch = []

        def process(rows):
            regions = []
            for number, fields in rows:
                if len(fields) < 12:
                    continue
                contig, start, end, strand = fields[5], int(fields[7]), int(fields[8]), fields[4]
                if contig not in lengths:
                    continue
                if not 0 <= start < end <= lengths[contig]:
                    continue
                regions.append('{}:{}-{}'.format(contig, start+1, end))
                regions.append('{}:{}-{}'.format(contig, end+1, end+3) if strand == '+'
                               else '{}:{}-{}'.format(contig, start-2, start))
            fetched = independent.fetch_regions(a.samtools, a.reference, regions) if regions else {}
            for number, fields in rows:
                counts['rows'] += 1
                target = pam = ''
                if len(fields) >= 12 and fields[5] in lengths:
                    start, end = int(fields[7]), int(fields[8])
                    if 0 <= start < end <= lengths[fields[5]]:
                        target = fetched.get('{}:{}-{}'.format(fields[5], start+1, end), '')
                        pam = (fetched.get('{}:{}-{}'.format(fields[5], end+1, end+3), '') if fields[4] == '+'
                               else fetched.get('{}:{}-{}'.format(fields[5], start-2, start), ''))
                        if fields[4] == '-' and target and pam:
                            target, pam = independent.rc(target), independent.rc(pam)
                result, accounted = label_row(number, fields, guides, target, pam, ops_cfd if target and pam else None)
                counts['accounted'] += accounted
                key = (fields[0], fields[5], fields[4], fields[8] if fields[4] == '+' else fields[7]) if len(fields) >= 12 else None
                result['contained_in_new'] = key in new_keys if key else False
                counts['rows_contained_in_new'] += result['contained_in_new']
                if result['strict_eligible']:
                    counts['strict_eligible'] += 1
                    counts['strict_eligible_contained_in_new'] += result['contained_in_new']
                    counts['strict_eligible_missing_in_new'] += not result['contained_in_new']
                if result['reconstructed_cfd'] != '':
                    counts['cfd_compared'] += 1
                    counts['cfd_disagreements'] += not result['cfd_match']
                writer.writerow(result)
        for number, line in enumerate(source, 1):
            batch.append((number, line.rstrip('\n').split('\t')))
            if len(batch) == 500:
                process(batch)
                batch = []
        if batch:
            process(batch)
    require_total = counts['rows']
    summary = dict(counts, status='PASS' if counts['accounted'] == counts['rows'] else 'FAIL',
                   guides=len(guides), new_anchor_keys=len(new_keys),
                   source_sha256=digest(a.paf), new_paf_sha256=digest(a.new_paf),
                   caveat='Legacy-row labeling only: strict eligible CIGAR witnesses are separate from '
                          'standards-valid PAF geometry and from genome-wide recall; old scores are not '
                          'reinterpreted and old rows are never repaired.')
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (out/'INCOMPLETE').unlink()
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
