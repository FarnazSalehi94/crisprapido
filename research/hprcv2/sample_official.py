#!/usr/bin/env python3
"""Prepare/check small independent official-inference samples; no search/training."""
import argparse
import csv
import itertools
import json
from pathlib import Path

from collect import collect, only
from inventory import digest


def prepare(root, project, out):
    report = collect(root, project)
    if report['status'] != 'COMPLETE_SAFETY_PILOT':
        raise ValueError('pilot must be complete')
    out.mkdir()
    samples = []
    for task, result in report['configurations'].items():
        marker = json.loads(Path(result['completion_marker']).read_text())
        groups, tails = {}, []
        with only(marker['outputs'], 'records.tsv').open() as f, only(marker['outputs'], 'annotations.tsv').open() as a:
            for audited, annotation in itertools.zip_longest(csv.DictReader(f, delimiter='\t'), csv.DictReader(a, delimiter='\t')):
                if audited is None or annotation is None or audited['line'] != str(int(annotation['row_index'])+1):
                    raise ValueError('audit/annotation identity mismatch')
                if annotation['crispr_bulge_reason'] != 'scored':
                    continue  # No official prediction outside supported policy.
                if any(audited[k] != annotation[k] for k in ['aligned_guide', 'aligned_target', 'pam']):
                    raise ValueError('independent aligned input differs')
                gap = 'DNA' if '-' in audited['aligned_guide'] else 'RNA' if '-' in audited['aligned_target'] else 'none'
                sample = dict(task=task, row_index=int(annotation['row_index']), guide_id=audited['guide'],
                              cigar=annotation['cigar'], guide=audited['aligned_guide']+'NGG',
                              target=audited['aligned_target']+audited['pam'], expected=float(annotation['crispr_bulge']),
                              stratum=[audited['guide'], audited['mismatches'], gap])
                groups.setdefault(tuple(sample['stratum']), sample)
                tails.append(sample)
        selected = {(r['task'], r['row_index']): r for r in groups.values()}
        tails.sort(key=lambda r: (r['expected'], r['row_index']))
        for row in tails[:4]+tails[-4:]:
            selected[(row['task'], row['row_index'])] = row
        samples.extend(selected[k] for k in sorted(selected))
    for n, row in enumerate(samples):
        row['id'] = n
    with (out/'requests.jsonl').open('x') as f:
        for row in samples:
            f.write(json.dumps(row)+'\n')
    with (out/'sample-provenance.json').open('x') as f:
        json.dump(dict(project=project, sample_rows=len(samples), requests_sha256=digest(out/'requests.jsonl'),
                       completions={k: v['completion_sha256'] for k, v in report['configurations'].items()},
                       sampling='First supported row per guide/mismatch/gap-type stratum plus four low/high prediction tails per configuration; not biological accuracy/calibration.',
                       input_reconstruction='Independent samtools/CIGAR audit, not adapter feature cropping or realignment.'), f, indent=2)
        f.write('\n')
    return len(samples)


def check(out):
    provenance = json.loads((out/'sample-provenance.json').read_text())
    if digest(out/'requests.jsonl') != provenance['requests_sha256']:
        raise ValueError('sample input changed')
    requested = [json.loads(s) for s in (out/'requests.jsonl').read_text().splitlines()]
    predicted = [json.loads(s) for s in (out/'official.jsonl').read_text().splitlines()]
    if [r['id'] for r in requested] != [r['id'] for r in predicted]:
        raise ValueError('official sample row count/identity mismatch')
    monitor = json.loads((out/'monitor-summary.json').read_text())
    if monitor['reason'] != 'completed':
        raise ValueError('official reference supervisor failed')
    errors = [abs(r['expected']-p['score']) for r, p in zip(requested, predicted)]
    if not errors or max(errors) > 1e-6 or any(not 0 <= p['score'] <= 1 for p in predicted):
        raise ValueError('official prediction disagreement')
    result = dict(status='PASS', rows=len(errors), max_absolute_difference=max(errors), tolerance=1e-6,
                  scope='Sampled implementation agreement, not biological accuracy, calibrated probability or genome-wide recall.',
                  official_metrics=json.loads((out/'official-metrics.json').read_text()), monitor=monitor,
                  input_sha256=digest(out/'requests.jsonl'), official_output_sha256=digest(out/'official.jsonl'))
    with (out/'agreement.json').open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', choices=['sbds', 'brunello'])
    p.add_argument('--root', type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    if a.check:
        print(json.dumps(check(a.out)))
    else:
        if not a.root or not a.project or a.root.name != a.project+'_pam_aware_v1':
            raise ValueError('project/root required')
        print(json.dumps(dict(samples=prepare(a.root, a.project, a.out))))


if __name__ == '__main__':
    main()
