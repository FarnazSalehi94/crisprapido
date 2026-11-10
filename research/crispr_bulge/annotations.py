"""Human-readable per-hit annotations. N/A is never a numeric zero."""
import re

from adapter import parse_paf, rc

COLUMNS = ['row_index', 'guide_id', 'contig', 'start_0based', 'end_0based', 'strand',
           'guide_sequence', 'target_sequence', 'pam', 'cigar', 'aligned_guide',
           'aligned_target', 'mismatches', 'bulge_events', 'max_bulge_size',
           'cfd', 'cfd_reason', 'crispr_bulge', 'crispr_bulge_reason']


def annotation(n, raw, guides, reference, state, score):
    f, tags = parse_paf(raw)
    guide = guides[f[0]]
    lo, hi = int(f[7]), int(f[8])
    target = reference.fetch(f[5], lo, hi)
    if f[4] == '-':
        target = rc(target)
    pam = 'N/A'
    if f[4] == '+' and hi+3 <= int(f[6]):
        pam = reference.fetch(f[5], hi, hi+3)
    elif f[4] == '-' and lo >= 3:
        pam = rc(reference.fetch(f[5], lo-3, lo))
    cigar = tags.get('cg', ('', 'N/A'))[1]
    ag, at = [], []
    qi, ti = int(f[2]), 0
    if not re.fullmatch(r'(?:[1-9][0-9]*[=XMID])+', cigar):
        raise ValueError('cannot annotate malformed CIGAR')
    for count, op in re.findall(r'(\d+)([=XMID])', cigar):
        for _ in range(int(count)):
            ag.append('-' if op == 'D' else guide[qi])
            at.append('-' if op == 'I' else target[ti])
            qi += op != 'D'; ti += op != 'I'
    if qi != int(f[3]) or ti != len(target):
        raise ValueError('annotation consumption mismatch')
    cfd = tags.get('cf', ('', 'N/A'))[1]
    return [n, f[0], f[5], lo, hi, f[4], guide, target, pam, cigar, ''.join(ag), ''.join(at),
            *[tags.get(k, ('', 'N/A'))[1] for k in ('nm', 'ng', 'bs')],
            cfd, 'available' if cfd != 'N/A' else 'upstream_cfd_unavailable',
            'N/A' if score is None else format(score, '.9g'), state]
