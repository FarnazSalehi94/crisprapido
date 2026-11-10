#!/usr/bin/env python3
"""Exercise all CLI modes on saved job 2871283, with full official score parity.

No genome search; never modify original input/reference/index/model files.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from adapter import digest, encode_tokens
from evaluate_assembly import audit, independent_requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def run(command, **kwargs):
    subprocess.run(list(map(str, command)), check=True, timeout=900, **kwargs)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('cli', 'input', 'reference', 'reference-index', 'upstream', 'model-python', 'out'):
        p.add_argument('--'+name, required=True, type=Path)
    a = p.parse_args(); a.out.mkdir(exist_ok=False, parents=True)
    original = digest(a.input)
    if original != 'a28d0574a915a57d431851ab681b0602b8f5fb1740d1c8efeb6e3b4a378f23be':
        raise ValueError('not saved job 2871283 output')
    reference_sha = 'bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b'
    guide = 'GCCAGGTACAGGTTCCGTCC'
    guides = a.out/'guides.fa'; guides.write_text('>pilot\n'+guide+'\n')
    results = {}; all_rows = {}
    for mode in ('cfd', 'crispr-bulge', 'both'):
        call = [a.cli, '-r', a.reference, '-G', guides, '--score-paf', a.input, '--scoring', mode]
        if mode != 'cfd':
            call += ['--scoring-output-dir', a.out/mode, '--crispr-bulge-python', a.model_python,
                     '--crispr-bulge-model-root', a.upstream, '--scoring-reference-sha256', reference_sha,
                     '--scoring-reference-index', a.reference_index, '--scoring-reference-index-sha256', digest(a.reference_index)]
        run([sys.executable, REPO/'research/pam_anchored_poc/monitored_pilot/supervise.py',
             '--limit-kib', '4194304', '--interval', '.02', '--log', a.out/(mode+'-tree.log'),
             '--summary', a.out/(mode+'-tree.json'), '--', sys.executable, HERE/'measure.py',
             '--metrics', a.out/(mode+'-measurement.json'), '--stdout', a.out/(mode+'.paf'),
             '--stderr', a.out/(mode+'.stderr'), '--timeout', '900', '--', *call], env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
        result = dict(measurement=json.loads((a.out/(mode+'-measurement.json')).read_text()),
                      tree=json.loads((a.out/(mode+'-tree.json')).read_text()), output_sha256=digest(a.out/(mode+'.paf')))
        if mode == 'cfd':
            if result['output_sha256'] != original:
                raise ValueError('default PAF changed')
            result['rows'] = len((a.out/(mode+'.paf')).read_text().splitlines())
        else:
            scored = a.out/mode/'scoring'
            rows, counts, reasons, geometry, _ = audit(a.input, scored/'output.paf')
            if counts != dict(records=21447, cfd=21447, crispr_bulge=4526, both=4526, crispr_bulge_unavailable=16921):
                raise ValueError('unexpected availability')
            if digest(scored/'output.paf') != result['output_sha256']:
                raise ValueError('CLI stdout/bundle mismatch')
            with (scored/'annotations.tsv').open() as f:
                annotations = list(csv.DictReader(f, delimiter='\t'))
            if len(annotations) != len(rows):
                raise ValueError('sidecar lost records')
            for r, annotation in zip(rows, annotations):
                f = r['fields']
                if (int(annotation['row_index']) != r['id'] or annotation['guide_sequence'] != guide
                    or annotation['contig'] != f[5] or annotation['start_0based'] != f[7]
                    or annotation['end_0based'] != f[8] or annotation['strand'] != f[4]
                    or annotation['cigar'] != r['cigar'] or float(annotation['cfd']) != r['cfd']
                    or annotation['crispr_bulge_reason'] != r['status']
                    or (r['cb'] is None and annotation['crispr_bulge'] != 'N/A')
                    or (r['cb'] is not None and float(annotation['crispr_bulge']) != r['cb'])):
                    raise ValueError('sidecar identity/score mismatch')
            result.update(counts=counts, unavailable_reasons=reasons, unavailable_geometry=geometry,
                          metadata=json.loads((scored/'metadata.json').read_text()), sidecar_rows=len(annotations))
            all_rows[mode] = rows
        results[mode] = result
    # Independently extract/replay EVERY supported representative, not just a sample.
    supported = [r for r in all_rows['both'] if r['cb'] is not None]
    work = a.out/'official-all'; work.mkdir()
    requests_path, requests, strata = independent_requests(supported, a.reference, work, guide)
    run([a.model_python, '-I', HERE/'official_reference.py', '--upstream', a.upstream,
         '--input', requests_path, '--output', work/'predictions.jsonl', '--metrics', work/'metrics.json'],
        stdout=(work/'stdout').open('x'), stderr=(work/'stderr').open('x'))
    official = [json.loads(s) for s in (work/'predictions.jsonl').read_text().splitlines()]
    if [r['id'] for r in official] != [r['id'] for r in supported]:
        raise ValueError('official IDs/count mismatch')
    import numpy as np
    tokens_sha = hashlib.sha256(np.asarray(encode_tokens(requests), dtype='<i8').tobytes()).hexdigest()
    metrics = json.loads((work/'metrics.json').read_text())
    if tokens_sha != metrics['categorical_sha256']:
        raise ValueError('official encoding mismatch')
    errors = {}
    for mode, rows in all_rows.items():
        scored = [r for r in rows if r['cb'] is not None]
        error = max(abs(r['cb']-o['score']) for r, o in zip(scored, official))
        if error > 1e-6:
            raise ValueError('official prediction disagreement')
        errors[mode] = error
    if digest(a.input) != original or digest(a.reference) != reference_sha:
        raise ValueError('original data changed')
    sources = [HERE/n for n in ('score_paf.py', 'worker.py', 'adapter.py', 'annotations.py', 'reference_index.py', 'validate_cli.py')]
    sources += [REPO/'src'/n for n in ('main.rs', 'scoring.rs', 'reporting.rs', 'cfd_score.rs')]
    summary = dict(baseline='90d958a', input_sha256=original, reference_sha256=reference_sha,
        index_sha256=digest(a.reference_index), original_data_unchanged=True, modes=results,
        official_validation=dict(rows=len(official), max_absolute_errors=errors, tolerance=1e-6,
            metrics=metrics, strata=strata, requests_sha256=digest(requests_path), predictions_sha256=digest(work/'predictions.jsonl')),
        cli_sha256=digest(a.cli), source_sha256={str(s.relative_to(REPO)): digest(s) for s in sources},
        scope='saved PAF scoring only; no search, model-domain extension, biological accuracy or recall claim')
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
