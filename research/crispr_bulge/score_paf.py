#!/usr/bin/env python3
"""Isolated PAF postprocessor: cfd (byte-preserving), crispr-bulge, or both.

CFD is calculated upstream by unchanged CRISPRapido, not reimplemented here.
This tool never changes representatives, search options, or --min-cfd behavior.
"""
import argparse
import csv
from collections import Counter
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adapter import IndexedFasta, convert, digest, guides_from_fasta, parse_paf
from annotations import COLUMNS, annotation
from reference_index import PinnedIndexedFasta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--out-dir', required=True, type=Path, help='NEW exclusive result directory')
    p.add_argument('--scoring', choices=('cfd', 'crispr-bulge', 'both'), default='cfd')
    p.add_argument('--reference', type=Path, help='plain, regularly wrapped FASTA')
    g = p.add_mutually_exclusive_group()
    g.add_argument('--guides', type=Path)
    g.add_argument('--guide', help='single-guide sequence; matches the legacy PAF query name Guide')
    p.add_argument('--expected-reference-sha256', help='optional trusted digest from the search run')
    p.add_argument('--model-root', type=Path)
    p.add_argument('--model-python', type=Path, help='absolute dedicated TensorFlow interpreter')
    p.add_argument('--reference-index', type=Path, help='explicit trusted prepared offset index')
    p.add_argument('--reference-index-sha256', help='required pin for prepared index')
    p.add_argument('--annotations', action='store_true', help='write annotations.tsv including N/A and sequence context')
    p.add_argument('--preserve-cfd', action='store_true', help='retain upstream CFD tags even in CB-only mode')
    p.add_argument('--batch-size', type=int, default=512)
    p.add_argument('--worker-timeout', type=float, default=300)
    a = p.parse_args()
    if not 1 <= a.batch_size <= 4096 or not math.isfinite(a.worker_timeout) or a.worker_timeout <= 0:
        p.error('invalid batch size/worker timeout')
    if bool(a.reference_index) != bool(a.reference_index_sha256):
        p.error('reference index and its SHA256 must be supplied together')
    if a.reference_index and not a.expected_reference_sha256:
        p.error('prepared index requires trusted expected reference SHA256')
    if a.annotations and a.scoring == 'cfd':
        p.error('annotations are only available in opted-in neural modes')
    if a.scoring != 'cfd':
        if not a.reference or not (a.guides or a.guide) or not a.model_root or not a.model_python:
            p.error('CRISPR-Bulge modes require reference, guide(s), model root and model Python')
        if not a.model_python.is_absolute():
            p.error('model Python must be an absolute path')
        if not a.model_root.is_dir():
            p.error('model root is missing')
    a.out_dir.mkdir(exist_ok=False, parents=True)
    status = a.out_dir/'status'
    status.write_text('INCOMPLETE\n')
    started = time.perf_counter()
    source_sha = digest(a.input)
    metadata = dict(scoring=a.scoring, input_sha256=source_sha,
                    min_cfd='not applied here; upstream filter, if any, is preserved',
                    representative_policy='input PAF representative unchanged',
                    rows=0, statuses={}, worker=None)
    output = a.out_dir/'output.paf.incomplete'
    reference = None
    try:
        if a.scoring == 'cfd':
            # Exact pass-through of upstream CFD values, formatting and ordering.
            # Reject accidental double annotation rather than labeling CB output CFD-only.
            with a.input.open(encoding='ascii', newline='') as source:
                for raw in source:
                    parse_paf(raw)
                    metadata['rows'] += 1
            shutil.copyfile(a.input, output)
            metadata['statuses'] = {'upstream_cfd_preserved': metadata['rows']}
        else:
            indexed = time.perf_counter()
            reference = (PinnedIndexedFasta(a.reference, a.reference_index, a.reference_index_sha256,
                                             a.expected_reference_sha256) if a.reference_index else IndexedFasta(a.reference))
            metadata['reference_index_sha256'] = a.reference_index_sha256
            metadata['reference_index_strategy'] = 'pinned_offsets_full_pre_post_reference_hash' if a.reference_index else 'validate_wrapping_alphabet_and_hash'
            if a.expected_reference_sha256 and reference.sha256 != a.expected_reference_sha256:
                raise ValueError('reference digest does not match expected search reference')
            metadata['reference_sha256'] = reference.sha256
            metadata['reference_index_scan_seconds'] = time.perf_counter()-indexed
            guide_sha = digest(a.guides) if a.guides else None
            guides = guides_from_fasta(a.guides) if a.guides else {'Guide': a.guide.upper()}
            metadata['single_guide'] = a.guide.upper() if a.guide else None
            metadata['model_root'] = str(a.model_root.resolve())
            request = a.out_dir/'requests.jsonl.incomplete'
            spool = a.out_dir/'records.jsonl.incomplete'
            counts = Counter()
            seen = set()
            conversion = time.perf_counter()
            with a.input.open(encoding='ascii', newline='') as source, request.open('x') as requests, spool.open('x') as records:
                for n, raw in enumerate(source):
                    key, state, features = convert(raw, guides, reference)
                    if key in seen:
                        raise ValueError('duplicate input PAM anchor')
                    seen.add(key)
                    counts[state] += 1
                    details = annotation(n, raw, guides, reference, state, None) if a.annotations else None
                    records.write(json.dumps(dict(id=n, raw=raw, status=state, annotation=details))+'\n')
                    if features is not None:
                        requests.write(json.dumps(dict(id=n, **features))+'\n')
            metadata['rows'] = sum(counts.values())
            metadata['statuses'] = dict(sorted(counts.items()))
            metadata['conversion_seconds'] = time.perf_counter()-conversion
            scores = a.out_dir/'scores.jsonl.incomplete'
            metrics = a.out_dir/'worker-metrics.json'
            manifest = Path(__file__).with_name('model_manifest.json')
            metadata['model_manifest_sha256'] = digest(manifest)
            metadata['model_id'] = json.loads(manifest.read_text())['model_id']
            if counts['scored']:
                call = [str(a.model_python), '-I', str(Path(__file__).with_name('worker.py').resolve()),
                        '--model-root', str(a.model_root.resolve()), '--input', str(request.resolve()),
                        '--output', str(scores.resolve()), '--metrics', str(metrics.resolve()),
                        '--batch-size', str(a.batch_size)]
                with (a.out_dir/'worker.stdout').open('x') as stdout, (a.out_dir/'worker.stderr').open('x') as stderr:
                    run = subprocess.run(call, stdout=stdout, stderr=stderr, timeout=a.worker_timeout,
                                         env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
                if run.returncode:
                    raise RuntimeError('model worker failed; see worker.stderr')
                metadata['worker'] = json.loads(metrics.read_text())
                if metadata['worker']['model_id'] != metadata['model_id']:
                    raise RuntimeError('worker model identity mismatch')
            else:
                scores.touch(exist_ok=False)
                metadata['model_not_loaded_reason'] = 'no inputs in conservative model domain'
            annotations = (a.out_dir/'annotations.tsv.incomplete').open('x', newline='') if a.annotations else None
            try:
                if annotations:
                    writer = csv.writer(annotations, delimiter='\t', lineterminator='\n')
                    writer.writerow(COLUMNS)
                write_results(scores, spool, output, a, annotations, writer if annotations else None)
            finally:
                if annotations:
                    annotations.close()
            if digest(manifest) != metadata['model_manifest_sha256']:
                raise ValueError('model manifest changed during scoring')
            if digest(a.reference) != reference.sha256 or (a.guides and digest(a.guides) != guide_sha):
                raise ValueError('reference/guides changed during scoring')
            if a.reference_index and digest(a.reference_index) != a.reference_index_sha256:
                raise ValueError('reference index changed during scoring')
            metadata['guide_sha256'] = guide_sha
        if digest(a.input) != source_sha:
            raise ValueError('input PAF changed during scoring')
        metadata['output_sha256'] = digest(output)
        metadata['elapsed_seconds'] = time.perf_counter()-started
        metadata['scope'] = 'post-search scoring only; CFD unchanged upstream; no genomic recall claim'
        if a.annotations:
            metadata['annotations_sha256'] = digest(a.out_dir/'annotations.tsv.incomplete')
        (a.out_dir/'metadata.json').write_text(json.dumps(metadata, indent=2, sort_keys=True)+'\n')
        output.rename(a.out_dir/'output.paf')
        if a.annotations:
            (a.out_dir/'annotations.tsv.incomplete').rename(a.out_dir/'annotations.tsv')
        status.write_text('COMPLETE scoring only; no search or recall validation\n')
        # Requests/spool are transient, not a reusable verified reference cache.
        for name in ('requests.jsonl.incomplete', 'records.jsonl.incomplete', 'scores.jsonl.incomplete'):
            path = a.out_dir/name
            if path.exists():
                path.unlink()
    finally:
        if reference:
            reference.close()


def write_results(scores, spool, output, a, annotations, writer):
    with scores.open() as predictions, spool.open() as records, output.open('x', encoding='ascii', newline='') as target:
        for line in records:
            row = json.loads(line)
            f, _ = parse_paf(row['raw'])
            if a.scoring == 'crispr-bulge' and not a.preserve_cfd:
                f = f[:12]+[text for text in f[12:] if not text.startswith('cf:')]
            score = None
            if row['status'] == 'scored':
                item = predictions.readline()
                if not item:
                    raise RuntimeError('model omitted an eligible record')
                item = json.loads(item)
                if item['id'] != row['id'] or type(item['score']) not in (float, int) or not math.isfinite(item['score']) or not 0 <= item['score'] <= 1:
                    raise RuntimeError('model output ID/value mismatch')
                score = item['score']
                f.append('cb:f:'+format(score, '.9g'))
            if annotations:
                details = row['annotation']
                details[-2] = 'N/A' if score is None else format(score, '.9g')
                writer.writerow(details)
            f.append('bu:Z:'+row['status'])
            eol = '\r\n' if row['raw'].endswith('\r\n') else '\n' if row['raw'].endswith('\n') else ''
            target.write('\t'.join(f)+eol)
        if predictions.read():
            raise RuntimeError('unexpected extra model output')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print('isolated CRISPR-Bulge scorer: '+str(e), file=sys.stderr)
        sys.exit(1)
