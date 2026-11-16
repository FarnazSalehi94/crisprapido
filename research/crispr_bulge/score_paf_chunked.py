#!/usr/bin/env python3
"""Memory-safe chunked wrapper around the validated score_paf.py scorer.

Drop-in replacement for the scoring script the CRISPRapido CLI invokes: same
argument interface, same out-dir contract (status/output.paf/annotations.tsv/
metadata.json). Splits the input PAF into bounded row chunks, scores each
chunk with the UNMODIFIED pinned score_paf.py (no model, CFD formula,
representative, or eligibility change), then merges with strict row-count,
ordering, and checksum verification.

Memory: the inner scorer's dominant consumer is its per-input duplicate-key
set, which scales with rows; chunking bounds it to --chunk-rows. Restart-safe:
verified completed chunks are skipped; chunk outputs are removed only after
the fully verified merge is published.

Global integrity: the inner scorer's duplicate-anchor check is chunk-local.
The merge instead enforces a strictly stronger end-to-end property: every
input row appears exactly once, in order, as the prefix of its annotated
output row (streaming, O(1) memory).
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import uuid

COLUMNS = 'row_index guide_id contig start_0based end_0based strand guide_sequence target_sequence pam cigar aligned_guide aligned_target mismatches bulge_events max_bulge_size cfd cfd_reason crispr_bulge crispr_bulge_reason'.split()


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def count_lines(path):
    n = 0
    with open(path, 'rb') as f:
        for _ in f:
            n += 1
    return n


def split_input(input_paf, chunks_dir, chunk_rows, resume=False):
    chunks_dir.mkdir(parents=True, exist_ok=True)
    chunks = []
    with open(input_paf, 'rb') as src:
        index = 0
        while True:
            first = src.readline()
            if not first:
                break
            path = chunks_dir/('chunk-%05d.paf' % index)
            if resume and path.exists():
                path.unlink()  # derived input; no live reader can exist on resume
            rows = 1
            with path.open('xb') as out:
                out.write(first)
                while rows < chunk_rows:
                    line = src.readline()
                    if not line:
                        break
                    out.write(line)
                    rows += 1
            chunks.append(dict(index=index, rows=rows, path=str(path),
                               sha256=digest(path)))
            index += 1
    if not chunks:
        raise ValueError('empty input PAF')
    return chunks


def chunk_out_dir(root, index):
    return root/('chunks/%05d' % index)


def score_chunk(args, inner, chunk, out_dir):
    # The inner scorer creates out_dir itself (exclusive); the wrapper keeps
    # its invocation log beside it.
    command = [sys.executable, '-I', '-B', str(inner),
        '--input', chunk['path'], '--out-dir', out_dir, '--scoring', args.scoring,
        '--reference', args.reference, '--model-root', args.model_root,
        '--model-python', args.model_python, '--batch-size', str(args.batch_size),
        '--worker-timeout', str(args.worker_timeout)]
    if args.annotations:
        command.append('--annotations')
    if args.preserve_cfd:
        command.append('--preserve-cfd')
    if args.guides:
        command += ['--guides', args.guides]
    elif args.guide:
        command += ['--guide', args.guide]
    if args.expected_reference_sha256:
        command += ['--expected-reference-sha256', args.expected_reference_sha256]
    if args.reference_index:
        command += ['--reference-index', args.reference_index,
                    '--reference-index-sha256', args.reference_index_sha256]
    before = resource.getrusage(resource.RUSAGE_CHILDREN); start = time.monotonic()
    log = out_dir.parent/('invocation-%05d.stdout' % chunk['index'])
    with log.open('x') as stdout:
        result = subprocess.run(command, stdout=stdout, stderr=subprocess.STDOUT, check=True)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    status = (out_dir/'status').read_text()
    if not status.startswith('COMPLETE'):
        raise RuntimeError('inner scorer did not complete chunk '+str(chunk['index']))
    metadata = json.loads((out_dir/'metadata.json').read_text())
    if metadata['rows'] != chunk['rows']:
        raise ValueError('chunk row-count mismatch: '+str(chunk['index']))
    if digest(out_dir/'output.paf') != metadata['output_sha256']:
        raise ValueError('chunk output hash mismatch: '+str(chunk['index']))
    annotation_rows = None
    if args.annotations:
        annotation_rows = count_lines(out_dir/'annotations.tsv') - 1  # header
        if annotation_rows != chunk['rows']:
            raise ValueError('chunk annotation row mismatch: '+str(chunk['index']))
    return dict(elapsed_seconds=round(time.monotonic()-start, 1),
                child_user_seconds=round(after.ru_utime-before.ru_utime, 1),
                child_peak_rss_kib=after.ru_maxrss, output_sha256=metadata['output_sha256'],
                annotation_rows=annotation_rows, returncode=result.returncode)


def merge(args, chunks, records, out_dir):
    output = out_dir/'output.paf.incomplete'
    annotations = out_dir/'annotations.tsv.incomplete' if args.annotations else None
    written = 0
    with output.open('xb') as target:
        for chunk, record in zip(chunks, records):
            path = chunk_out_dir(out_dir, chunk['index'])/'output.paf'
            if digest(path) != record['output_sha256']:
                raise ValueError('chunk output changed before merge: '+str(chunk['index']))
            with path.open('rb') as f:
                shutil.copyfileobj(f, target, 8*1024*1024)
            written += chunk['rows']
    if annotations:
        global_index = 0
        with annotations.open('x', newline='') as target:
            target.write('\t'.join(COLUMNS)+'\n')
            for chunk in chunks:
                with (chunk_out_dir(out_dir, chunk['index'])/'annotations.tsv').open() as f:
                    f.readline()  # per-chunk header
                    for line in f:
                        fields = line.rstrip('\n').split('\t')
                        fields[0] = str(global_index)  # restore global row order
                        target.write('\t'.join(fields)+'\n')
                        global_index += 1
    return written


def verify_merged(input_paf, merged, expected_rows):
    # Every input row must appear exactly once, in order, as the prefix of
    # its annotated output row (one extra trailing tag field).
    n = 0
    with open(input_paf, 'rb') as src, open(merged, 'rb') as out:
        for raw in src:
            scored = out.readline()
            if not scored or not scored.startswith(raw.rstrip(b'\r\n') + b'\t'):
                raise ValueError('merged output diverges from input at row %d' % n)
            n += 1
        if out.read(1) != b'' or n != expected_rows:
            raise ValueError('merged output row count mismatch')
    return n


def publish(out_dir, metadata):
    for name in ('metadata.json', 'output.paf', 'annotations.tsv'):
        target = out_dir/name
        temporary = out_dir/(name+'.tmp-'+uuid.uuid4().hex)
        (out_dir/(name+'.incomplete')).replace(temporary)
        os.replace(temporary, target)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--out-dir', required=True, type=Path, help='NEW exclusive result directory')
    p.add_argument('--scoring', choices=('cfd', 'crispr-bulge', 'both'), default='cfd')
    p.add_argument('--reference', type=Path)
    p.add_argument('--guides', type=Path)
    p.add_argument('--guide')
    p.add_argument('--expected-reference-sha256')
    p.add_argument('--model-root', type=Path)
    p.add_argument('--model-python', type=Path)
    p.add_argument('--reference-index', type=Path)
    p.add_argument('--reference-index-sha256')
    p.add_argument('--annotations', action='store_true')
    p.add_argument('--preserve-cfd', action='store_true')
    p.add_argument('--batch-size', type=int, default=512)
    p.add_argument('--worker-timeout', type=float, default=300)
    p.add_argument('--chunk-rows', type=int, default=5_000_000,
                   help='rows per chunk; bounds the inner scorer memory')
    p.add_argument('--inner-script', type=Path,
                   default=Path(__file__).with_name('score_paf.py'),
                   help='the UNMODIFIED validated scorer this wrapper calls per chunk')
    a = p.parse_args()
    if not 1 <= a.chunk_rows <= 100_000_000:
        p.error('--chunk-rows must be 1..100000000')
    if a.scoring != 'cfd' and not (a.reference and (a.guides or a.guide)
                                  and a.model_root and a.model_python):
        p.error('CB modes require reference, guide(s), model root and model Python')
    started = time.monotonic()
    resume = False
    if a.out_dir.exists():
        status_path = a.out_dir/'status'
        resumable = (status_path.is_file() and status_path.read_text().startswith('INCOMPLETE')
                     and (a.out_dir/'chunks').is_dir())
        if not resumable:
            p.error('out-dir exists and is not a resumable INCOMPLETE chunked run')
        resume = True
    else:
        a.out_dir.mkdir(exist_ok=False)
    status_path = a.out_dir/'status'
    if not resume:
        status_path.write_text('INCOMPLETE chunked\n')
    metadata = dict(scoring=a.scoring, chunked=True, chunk_rows=a.chunk_rows,
        inner_script=str(a.inner_script), inner_script_sha256=digest(a.inner_script),
        input_sha256=digest(a.input), rows=0, chunks=[], scope=\
        'chunked post-search scoring; identical formulas/model/eligibility to the '
        'validated single-pass scorer; merged output verified row-for-row against input')
    chunks_dir = a.out_dir/'chunks'
    chunks_dir.mkdir(exist_ok=resume)
    chunks = split_input(a.input, chunks_dir, a.chunk_rows, resume=resume)
    metadata['rows'] = sum(c['rows'] for c in chunks)
    for chunk in chunks:
        out = chunk_out_dir(a.out_dir, chunk['index'])
        record = None
        if (out/'status').exists() and (out/'metadata.json').is_file():
            try:
                prior = json.loads((out/'metadata.json').read_text())
                if prior.get('rows') == chunk['rows'] and digest(out/'output.paf') == prior['output_sha256']:
                    record = dict(reused=True, elapsed_seconds=None, child_peak_rss_kib=None,
                                  output_sha256=prior['output_sha256'],
                                  annotation_rows=count_lines(out/'annotations.tsv')-1 if a.annotations else None)
            except (OSError, KeyError, ValueError):
                record = None  # unverifiable partial attempt: rescore chunk
        if record is None:
            if out.exists():
                shutil.rmtree(out)  # only this chunk's uncommitted scratch
            record = score_chunk(a, a.inner_script, chunk, out)
            record['reused'] = False
        record.update(index=chunk['index'], rows=chunk['rows'],
                      input_sha256=chunk['sha256'])
        metadata['chunks'].append(record)
        (a.out_dir/'metadata.json.incomplete').write_text(json.dumps(metadata, indent=2, sort_keys=True)+'\n')
    written = merge(a, chunks, metadata['chunks'], a.out_dir)
    verified_rows = verify_merged(a.input, a.out_dir/'output.paf.incomplete', metadata['rows'])
    if verified_rows != written:
        raise ValueError('merge/verify row disagreement')
    metadata['output_sha256'] = digest(a.out_dir/'output.paf.incomplete')
    if a.annotations:
        metadata['annotations_rows'] = count_lines(a.out_dir/'annotations.tsv.incomplete') - 1
        metadata['annotations_sha256'] = digest(a.out_dir/'annotations.tsv.incomplete')
        if metadata['annotations_rows'] != metadata['rows']:
            raise ValueError('merged annotation row mismatch')
    metadata['elapsed_seconds'] = round(time.monotonic()-started, 1)
    (a.out_dir/'metadata.json.incomplete').write_text(json.dumps(metadata, indent=2, sort_keys=True)+'\n')
    publish(a.out_dir, metadata)
    status_path.write_text('COMPLETE chunked; verified merge identical to input rows\n')
    # Transient chunk working data: derived deterministically from the input.
    # Removed only after the verified merge is published and checkpointed.
    for chunk in chunks:
        Path(chunk['path']).unlink(missing_ok=True)
        shutil.rmtree(chunk_out_dir(a.out_dir, chunk['index']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
