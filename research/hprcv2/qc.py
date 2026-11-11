"""Pilot accounting; independent sequence/CFD audit is required before this step."""
import collections
import csv
import itertools
import json
import math
from pathlib import Path
import zlib


def require(condition, message):
    if not condition:
        raise ValueError(message)


def tags(line):
    result = {}
    for item in line.rstrip('\n').split('\t')[12:]:
        key, kind, value = item.split(':', 2)
        require(key not in result, 'duplicate tag')
        result[key] = (kind, value)
    return result


def summarize(raw, annotated, annotations, audited, guide_ids):
    per_guide = {name: collections.Counter() for name in guide_ids}
    statuses, seen = collections.Counter(), set()
    compressors = {key: zlib.compressobj(6, zlib.DEFLATED, 31) for key in ('none', '0.05', '0.20')}
    compressed = collections.Counter()
    with open(raw) as a, open(annotated) as b, open(annotations) as c, open(audited) as d:
        annotation_rows = csv.DictReader(c, delimiter='\t')
        audit_rows = csv.DictReader(d, delimiter='\t')
        for n, group in enumerate(itertools.zip_longest(a, b, annotation_rows, audit_rows), 1):
            old, new, ann, audit = group
            require(all(x is not None for x in group), 'row loss/addition or annotation/audit mismatch')
            f, g = old.rstrip('\n').split('\t'), new.rstrip('\n').split('\t')
            require(g[:len(f)] == f, 'original PAF fields/tags/order changed')
            require(f[0] in per_guide and int(audit['line']) == n, 'unknown guide or audit order')
            t = tags(new)
            require(t.get('bu', ('', ''))[0] == 'Z', 'missing neural status')
            status = t['bu'][1]
            require(('cb' in t) == (status == 'scored'), 'neural availability/status inconsistency')
            require(ann['guide_id'] == f[0] and ann['contig'] == f[5] and ann['strand'] == f[4]
                    and ann['start_0based'] == f[7] and ann['end_0based'] == f[8]
                    and ann['cigar'] == t['cg'][1] and int(ann['row_index']) == n-1, 'annotation identity mismatch')
            for field in ['pam', 'target_sequence', 'aligned_guide', 'aligned_target']:
                require(ann[field] == audit['target' if field == 'target_sequence' else field], 'independent sequence/annotation mismatch: '+field)
            require(ann['crispr_bulge_reason'] == status, 'annotation unavailable reason mismatch')
            if 'cb' in t:
                score = float(t['cb'][1])
                require(t['cb'][0] == 'f' and math.isfinite(score) and 0 <= score <= 1, 'invalid neural score')
                require(float(ann['crispr_bulge']) == score, 'annotation prediction mismatch')
            else:
                require(ann['crispr_bulge'] == 'N/A' and status != 'scored', 'unavailable prediction fabricated')
            cf = t.get('cf')
            require(ann['cfd'] == ('N/A' if cf is None else cf[1]), 'annotation CFD mismatch')
            require(ann['cfd_reason'] == ('upstream_cfd_unavailable' if cf is None else 'available'), 'CFD reason mismatch')
            key = (f[0], f[5], f[4], f[8] if f[4] == '+' else f[7])
            require(key not in seen, 'duplicate guide/contig/strand/PAM-anchor')
            seen.add(key)
            counter = per_guide[f[0]]
            counter['eligible_hits'] += 1
            counter['cfd_defined'] += cf is not None
            counter['cb_defined'] += 'cb' in t
            counter['cb_unavailable'] += 'cb' not in t
            counter['perfect_20eq'] += t['cg'][1] == '20='
            counter['mismatches_'+audit['mismatches']] += 1
            counter['bulge_events_'+audit['bulges']] += 1
            counter['strand_'+f[4]] += 1
            statuses[status] += 1
            # Saved PAF CFDs are rounded; prioritization uses audited full-precision
            # reconstruction, NOT an implicit search filter or neural cutoff.
            value = None if audit['independent_cfd'] in ('', 'None') else float(audit['independent_cfd'])
            require((value is None) == (cf is None), 'CFD availability differs from independent audit')
            for label, threshold in [('none', None), ('0.05', .05), ('0.20', .20)]:
                keep = threshold is None or value is None or value >= threshold
                counter['view_'+label+'_defined_above'] += value is not None and (threshold is None or value >= threshold)
                counter['view_'+label+'_unscored_preserved'] += value is None
                counter['view_'+label+'_rows'] += keep
                if keep:
                    counter['view_'+label+'_uncompressed_paf_bytes'] += len(new.encode())
                    compressed[label] += len(compressors[label].compress(new.encode()))
    for key, compressor in compressors.items():
        compressed[key] += len(compressor.flush())
    return dict(status='PASS', per_guide={k: dict(v) for k, v in per_guide.items()},
                cb_statuses=dict(statuses), eligible_hits=sum(c['eligible_hits'] for c in per_guide.values()),
                analytical_view_gzip_bytes=dict(compressed),
                views='Derived only; full eligible PAF is retained. Unavailable CFD remains in a separate preserved category. Neural unavailability never filters a row.',
                thresholds='Independent reconstructed full-precision CFD; not four-decimal stored PAF and not biological activity/calibration.',
                caveat='Emitted-row soundness, identity and score-availability checks; not genome-wide recall or biological prediction validation.')
