#!/usr/bin/env python3
"""Audit an existing both-mode result and independently check a stratified sample.

No search, input rewriting, model-domain extension or CLI integration. Large
artifacts stay in the supplied external work directory. Requires scipy/numpy
for statistics and samtools for independent indexed reference extraction.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import re
import subprocess
import sys

from adapter import digest, encode_tokens


def alignment(cigar):
    ops = [(int(n), op) for n, op in re.findall(r'(\d+)([=XMID])', cigar)]
    if ''.join(str(n)+op for n, op in ops) != cigar:
        raise ValueError('unsupported alignment representation')
    ins = sum(n for n, op in ops if op == 'I')
    dele = sum(n for n, op in ops if op == 'D')
    kind = 'mixed' if ins and dele else 'RNA' if ins else 'DNA' if dele else 'ungapped'
    return ops, kind, ins, dele


def audit(source, annotated):
    originals = source.read_text().splitlines(keepends=True)
    outputs = annotated.read_text().splitlines(keepends=True)
    if len(originals) != len(outputs):
        raise ValueError('row count changed')
    counts = Counter(); reasons = Counter(); reasons_geometry = defaultdict(Counter)
    geometries = Counter(); rows = []
    for i, (raw, out) in enumerate(zip(originals, outputs)):
        base = raw.rstrip('\r\n')
        if not out.startswith(base+'\t') or not out.endswith('\r\n' if raw.endswith('\r\n') else '\n'):
            raise ValueError('original fields/order changed')
        fields = base.split('\t')
        tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in fields[12:]}
        added = out.rstrip('\r\n')[len(base)+1:].split('\t')
        extra = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in added}
        if set(extra) not in ({'bu'}, {'bu', 'cb'}) or len(extra) != len(added):
            raise ValueError('unexpected annotations')
        status = extra['bu']; cb = float(extra['cb']) if 'cb' in extra else None
        cfd = float(tags['cf']) if 'cf' in tags else None
        if (status == 'scored') != (cb is not None):
            raise ValueError('score/status mismatch')
        for value in (cfd, cb):
            if value is not None and not (math.isfinite(value) and 0 <= value <= 1):
                raise ValueError('invalid score')
        ops, kind, ins, dele = alignment(tags['cg'])
        geometry = f'{kind}:I{ins}:D{dele}'
        # Independently cross-check first-failure domain reasons from CIGAR/nm.
        # Sequence/PAM ambiguity is an earlier check in the scoring adapter.
        if status in {'scored', 'alignment_exceeds_24_columns', 'unsupported_mixed_bulge',
                      'unsupported_multiple_bulge_events', 'unsupported_multibase_bulge',
                      'outside_documented_mismatch_domain'}:
            gap_events = sum(op in 'ID' and (j == 0 or ops[j-1][1] not in 'ID')
                             for j, (_, op) in enumerate(ops))
            expected = ('alignment_exceeds_24_columns' if sum(n for n, _ in ops)+3 > 24 else
                        'unsupported_mixed_bulge' if ins and dele else
                        'unsupported_multiple_bulge_events' if gap_events > 1 else
                        'unsupported_multibase_bulge' if ins+dele > 1 else
                        'outside_documented_mismatch_domain' if int(tags['nm']) > (4 if ins+dele else 6) else
                        'scored')
            if status != expected:
                raise ValueError('CIGAR/mismatch domain reason disagrees')
        geometries[geometry] += 1
        counts['records'] += 1
        counts['cfd'] += cfd is not None
        counts['crispr_bulge'] += cb is not None
        counts['both'] += cb is not None and cfd is not None
        counts['crispr_bulge_unavailable'] += cb is None
        if cb is None:
            reasons[status] += 1
            reasons_geometry[status][geometry] += 1
        rows.append(dict(id=i, fields=fields, cigar=tags['cg'], kind=kind,
                         mismatch=int(tags['nm']), cfd=cfd, cb=cb, status=status))
    return rows, dict(counts), dict(reasons), {k: dict(v) for k, v in reasons_geometry.items()}, dict(geometries)


def sample_rows(rows):
    """Four per observed kind/strand/mismatch stratum plus both score tails."""
    groups = defaultdict(list)
    supported = [r for r in rows if r['cb'] is not None]
    for row in supported:
        groups[(row['kind'], row['fields'][4], row['mismatch'])].append(row)
    rng = random.Random(2871283)
    chosen = {}
    for key in sorted(groups):
        for row in rng.sample(groups[key], min(4, len(groups[key]))):
            chosen[row['id']] = row
    for score in ('cfd', 'cb'):
        ranked = sorted((r for r in supported if r[score] is not None), key=lambda r: (r[score], r['id']))
        for row in ranked[:10]+ranked[-10:]:
            chosen[row['id']] = row
    return [chosen[i] for i in sorted(chosen)]


def compare_scores(rows):
    import numpy as np
    from scipy import stats
    paired = [r for r in rows if r['cfd'] is not None and r['cb'] is not None]
    def describe(group):
        a = np.array([[r['cfd'], r['cb']] for r in group])
        result = {'n': len(group)}
        if not len(group):
            return result
        for j, name in enumerate(('cfd', 'crispr_bulge')):
            result[name] = dict(zip(('min', 'p25', 'median', 'p75', 'p90', 'p99', 'max'),
                                   map(float, np.quantile(a[:, j], [0, .25, .5, .75, .9, .99, 1]))))
            result[name]['zeros'] = int((a[:, j] == 0).sum())
        if len(group) > 2 and all(np.ptp(a[:, j]) > 0 for j in (0, 1)):
            result.update(pearson=float(stats.pearsonr(a[:, 0], a[:, 1])[0]),
                          spearman=float(stats.spearmanr(a[:, 0], a[:, 1])[0]),
                          kendall_tau_b=float(stats.kendalltau(a[:, 0], a[:, 1])[0]))
        return result
    return {'paired': describe(paired), 'by_alignment': {
        k: describe([r for r in paired if r['kind'] == k]) for k in sorted({r['kind'] for r in paired})},
        'interpretation': 'Descriptive association on supported representatives only; CFD rounded to four decimals. No calibration, biological accuracy or interchangeability inference.'}


def independent_requests(sample, reference, work, guide):
    """samtools/htslib extraction + separately implemented CIGAR reconstruction.

    All faidx writes go to work, never beside the original reference.
    """
    index = work/'reference.fai'
    subprocess.run(['samtools', 'faidx', '--fai-idx', str(index), str(reference)], check=True)
    regions = work/'sample-regions.txt'
    spans = []
    for r in sample:
        f = r['fields']; lo, hi = int(f[7]), int(f[8])
        spans.append(f'{f[5]}:{lo+1 if f[4] == "+" else lo-2}-{hi+3 if f[4] == "+" else hi}')
    regions.write_text('\n'.join(spans)+'\n')
    extracted = subprocess.check_output(['samtools', 'faidx', '--fai-idx', str(index),
        '-r', str(regions), str(reference)], text=True)
    seqs = []
    for line in extracted.splitlines():
        if line.startswith('>'):
            seqs.append('')
        else:
            seqs[-1] += line.upper()
    if len(seqs) != len(sample):
        raise ValueError('independent extraction count mismatch')
    requests = []; strata = Counter()
    for r, seq in zip(sample, seqs):
        f = r['fields']
        if f[4] == '-':
            seq = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        target, pam = seq[:-3], seq[-3:]
        if set(seq)-set('ACGT') or pam[1:] != 'GG':
            raise ValueError('sample outside resolved NGG domain')
        qi = ti = mm = 0; ga = []; ta = []
        for n, op in alignment(r['cigar'])[0]:
            for _ in range(n):
                a = '-' if op == 'D' else guide[qi]
                b = '-' if op == 'I' else target[ti]
                qi += op != 'D'; ti += op != 'I'
                if op not in 'ID':
                    if (op == '=' and a != b) or (op == 'X' and a == b):
                        raise ValueError('independent CIGAR assertion failed')
                    mm += a != b
                ga.append(a); ta.append(b)
        if qi != len(guide) or ti != len(target) or mm != r['mismatch']:
            raise ValueError('independent consumption/mismatch check failed')
        requests.append(dict(id=r['id'], guide=''.join(ga)+'NGG', target=''.join(ta)+pam))
        strata[f'{r["kind"]}:{f[4]}:mm{mm}:{pam}'] += 1
    path = work/'sample-requests.jsonl'
    path.write_text(''.join(json.dumps(r)+'\n' for r in requests))
    return path, requests, dict(sorted(strata.items()))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'scored-dir', 'reference', 'work', 'upstream', 'model-python', 'output'):
        p.add_argument('--'+name, required=True, type=Path)
    p.add_argument('--guide', default='GCCAGGTACAGGTTCCGTCC')
    a = p.parse_args(); a.work.mkdir(exist_ok=False, parents=True)
    before = digest(a.input)
    if before != 'a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be':
        raise ValueError('not saved job 2871283 output')
    metadata = json.loads((a.scored_dir/'metadata.json').read_text())
    if metadata['input_sha256'] != before or metadata['output_sha256'] != digest(a.scored_dir/'output.paf'):
        raise ValueError('scoring artifact identity mismatch')
    rows, counts, reasons, reason_geometry, geometry = audit(a.input, a.scored_dir/'output.paf')
    sample = sample_rows(rows)
    requests_path, requests, strata = independent_requests(sample, a.reference, a.work, a.guide)
    here = Path(__file__).resolve().parent
    call = [str(a.model_python), '-I', str(here/'official_reference.py'), '--upstream', str(a.upstream),
            '--input', str(requests_path.resolve()), '--output', str((a.work/'official.jsonl').resolve()),
            '--metrics', str((a.work/'official-metrics.json').resolve())]
    subprocess.run([sys.executable, str(here/'measure.py'), '--metrics', str(a.work/'measurement.json'),
        '--stdout', str(a.work/'official.stdout'), '--stderr', str(a.work/'official.stderr'), '--', *call], check=True)
    official = [json.loads(s) for s in (a.work/'official.jsonl').read_text().splitlines()]
    if [r['id'] for r in official] != [r['id'] for r in sample]:
        raise ValueError('reference IDs/count changed')
    errors = [abs(r['cb']-o['score']) for r, o in zip(sample, official)]
    import numpy as np
    encoding_hash = hashlib.sha256(np.asarray(encode_tokens(requests), dtype='<i8').tobytes()).hexdigest()
    official_metrics = json.loads((a.work/'official-metrics.json').read_text())
    if encoding_hash != official_metrics['categorical_sha256'] or max(errors) > 1e-6:
        raise ValueError('official parity failed')
    if digest(a.input) != before:
        raise ValueError('original input changed')
    result = dict(input_sha256=before, original_unchanged=True, counts=counts,
        unavailable_reasons=reasons, unavailable_reason_by_geometry=reason_geometry,
        all_alignment_geometry=geometry, comparison=compare_scores(rows), scoring_metadata=metadata,
        sample_validation=dict(n=len(sample), selection='seed 2871283; up to four per kind/strand/mismatch stratum plus bottom/top ten of each score (ID tie-break)',
            strata=strata, max_absolute_error=max(errors), tolerance=1e-6,
            categorical_sha256=encoding_hash, official=official_metrics,
            request_sha256=digest(requests_path), predictions_sha256=digest(a.work/'official.jsonl'),
            reconstruction='independent samtools/htslib extraction and separate CIGAR replay'),
        software=dict(samtools=subprocess.check_output(['samtools', '--version'], text=True).splitlines()[0],
                      source_sha256={f.name: digest(f) for f in (here/'evaluate_assembly.py', here/'score_paf.py', here/'adapter.py', here/'worker.py', here/'official_reference.py')}),
        scope='one existing guide/assembly result; no search, biological accuracy or recall evaluation')
    with a.output.open('x') as f:
        json.dump(result, f, indent=2, sort_keys=True, allow_nan=False); f.write('\n')


if __name__ == '__main__':
    main()
