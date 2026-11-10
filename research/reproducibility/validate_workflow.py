#!/usr/bin/env python3
"""Bounded end-to-end PAM-anchored search/report/CFD/neural validation.

Uses existing fixture definitions, no genome input option. Never enables pruning
or score filtering. All indexes/output are created under a NEW external root.
"""
import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
BULGE = ROOT/'research/crispr_bulge'
sys.path[:0] = [str(BULGE), str(ROOT/'research/pam_anchored_poc/cfd_feasibility')]
from adapter import digest, guides_from_fasta
from evaluate_assembly import audit
from study import fixtures, load_tables, oracle, representatives, independent_cfd
from validate import independent_features, fasta_dict


def run(command, **kwargs):
    subprocess.run(list(map(str, command)), check=True, timeout=300, cwd=ROOT, **kwargs)


def parse_rows(text):
    """Reject duplicate anchors (query, contig, strand, PAM-adjacent boundary)."""
    rows = {}
    for raw in text.splitlines():
        f = raw.split('\t')
        tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in f[12:]}
        key = (f[0], f[5], f[4], int(f[8]) if f[4] == '+' else int(f[7]))
        if key in rows:
            raise ValueError('duplicate PAM anchor: '+str(key))
        rows[key] = (f, tags)
    return rows


def check_oracle(text, expected):
    actual = parse_rows(text)
    if set(actual) != set(expected):
        raise ValueError(f'anchor mismatch: missing={len(set(expected)-set(actual))}, unexpected={len(set(actual)-set(expected))}')
    for key, (f, tags) in actual.items():
        hit = expected[key]
        if (int(f[7]), int(f[8]), tags['cg'], int(tags['nm'])) != (hit.start, hit.end, hit.cigar, hit.mm):
            raise ValueError('representative mismatch')
        value = float(tags['cf']) if 'cf' in tags else None
        if (value is None) != (hit.score is None) or (value is not None and abs(value-hit.score) > .000051):
            raise ValueError('independent oracle CFD mismatch')
    return len(actual)


def independent_cfd_check(text, guides, refs, tables):
    counts = Counter()
    for raw in text.splitlines():
        features, _ = independent_features(raw, guides, refs)
        f = raw.split('\t'); tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in f[12:]}
        import re
        ops = [(int(n), op) for n, op in re.findall(r'(\d+)([=XMID])', tags['cg'])]
        expected = independent_cfd(features['guide'][:-3], features['target'][:-3], ops,
                                   features['target'][-3:], *tables)
        actual = float(tags['cf']) if 'cf' in tags else None
        if (actual is None) != (expected is None) or (actual is not None and abs(actual-expected) > .000051):
            raise ValueError('independent CFD replay mismatch')
        counts['defined' if actual is not None else 'unavailable'] += 1
    return dict(counts)


def sidecar_check(path, original, guides, refs):
    with path.open() as f:
        annotations = list(csv.DictReader(f, delimiter='\t'))
    rows = original.splitlines()
    if len(annotations) != len(rows):
        raise ValueError('annotation row loss')
    for n, (item, raw) in enumerate(zip(annotations, rows)):
        features, _ = independent_features(raw, guides, refs)
        f = raw.split('\t')
        if (int(item['row_index']) != n or item['guide_id'] != f[0]
                or item['guide_sequence'] != guides[f[0]]
                or item['aligned_guide'] != features['guide'][:-3]
                or item['aligned_target'] != features['target'][:-3]
                or item['target_sequence'] != features['target'][:-3].replace('-', '')
                or item['pam'] != features['target'][-3:]):
            raise ValueError('independent sidecar reconstruction mismatch')
        tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in f[12:]}
        if item['cfd'] != tags.get('cf', 'N/A'):
            raise ValueError('sidecar CFD changed')
        if (item['crispr_bulge'] == 'N/A') != (item['crispr_bulge_reason'] != 'scored'):
            raise ValueError('N/A/status mismatch')
    return annotations


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('cli', 'old-cli', 'native', 'builder', 'model-root', 'model-python', 'out'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--retained-panel', type=Path, help='optional pinned retained172-inputs directory')
    p.add_argument('--official-python', type=Path)
    p.add_argument('--official-root', type=Path)
    a = p.parse_args()
    if bool(a.official_python) != bool(a.official_root):
        p.error('both official reference options are required together')
    a.out.mkdir(parents=True, exist_ok=False)
    tables = load_tables(); guide, resolved, nrefs = fixtures()
    panels = [('resolved', resolved, {'g': guide}, 5, None),
              ('N', nrefs, {'g': guide}, 5, None),
              ('mixed_only', resolved, {'g': guide}, 0, None)]
    if a.retained_panel:
        expected_hashes = {'reference.fa': 'd64192d9b1cf365fc104bf8c9709f8566995fdddc64f3cac1930011f718df682',
            'guides.fa': 'ac5bc0038ac29ffa064c7fc704f3709bff08d5f10f4a40c1c6237baeab56a80c',
            'input.paf': 'a89f36435f4bfb388152740cff036430ab4bd3a824a20986c0a3701e0e53c943'}
        for name, sha in expected_hashes.items():
            if digest(a.retained_panel/name) != sha:
                raise ValueError('retained independent panel identity mismatch')
        panels.append(('retained172', fasta_dict(a.retained_panel/'reference.fa'),
                       guides_from_fasta(a.retained_panel/'guides.fa'), 5,
                       (a.retained_panel/'input.paf').read_text()))
    summary = dict(baseline='6fa96d3', host=platform.node(), panels=[], benchmark=[],
        scope='full bounded PAM-anchored CLI workflow; no whole-genome search; no filtering/pruning',
        python=platform.python_version(), official_check_requested=bool(a.official_python))
    sources = [a.cli, a.old_cli, a.native, a.builder, Path(__file__), BULGE/'adapter.py',
               BULGE/'worker.py', BULGE/'score_paf.py', ROOT/'research/pam_anchored_poc/poc.cpp',
               ROOT/'research/pam_anchored_poc/cfd_feasibility/study.py', ROOT/'mismatch_scores.txt', ROOT/'pam_scores.txt']
    pinned = {str(f): digest(f) for f in sources}
    official_requests = []; comparisons = []
    for label, refs, guides, mismatches, retained in panels:
        folder = a.out/label; folder.mkdir()
        ref = folder/'reference.fa'; guidefile = folder/'guides.fa'; idx = folder/'index'
        ref.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in refs.items()))
        guidefile.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in guides.items()))
        with (folder/'index-build.log').open('x') as log:
            run([sys.executable, ROOT/'research/pam_anchored_poc/compatibility/index_manifest.py',
                 'build', '--builder', a.builder, '--reference', ref, '--index', idx], stdout=log, stderr=subprocess.STDOUT)
        expected = None
        if retained is None:
            expected = {(name, *key): hit for name, g in guides.items()
                        for key, hit in representatives(oracle(refs, g, 'GG', mismatches, 2, .75, tables)).items()}
        command = [a.cli, '--search-backend', 'pam-anchored', '--pam-search-bin', a.native,
                   '--columba-index', idx, '-r', ref, '-G', guidefile, '-t', '1']
        # Exercise default mismatch/bulge/fraction/PAM settings, except explicit mixed-only control.
        if mismatches != 5:
            command += ['-m', str(mismatches)]
        old = subprocess.run(list(map(str, [a.old_cli, *command[1:]])), capture_output=True, text=True,
                             check=True, timeout=300, cwd=ROOT).stdout
        if retained is not None and parse_rows(old) != parse_rows(retained):
            raise ValueError('retained 172 representatives differ')
        if expected is not None:
            check_oracle(old, expected)
        cfd_counts = independent_cfd_check(old, guides, refs, tables)
        base = folder/'baseline.paf'; base.write_text(old)
        rows0 = parse_rows(old); statuses = None; results = []
        for repetition in range(3 if label == 'N' else 1):
            order = ['cfd', 'crispr-bulge', 'both']
            if repetition % 2:
                order.reverse()
            for mode in order:
                prefix = folder/f'{mode}-{repetition}'
                call = command + ['--scoring', mode]
                if mode != 'cfd':
                    call += ['--scoring-output-dir', prefix, '--crispr-bulge-python', a.model_python,
                             '--crispr-bulge-model-root', a.model_root,
                             '--scoring-reference-sha256', digest(ref)]
                with Path(str(prefix)+'.supervisor.stdout').open('x') as log:
                    run([sys.executable, ROOT/'research/pam_anchored_poc/monitored_pilot/supervise.py',
                        '--limit-kib', '4194304', '--interval', '.02', '--log', str(prefix)+'.rss',
                        '--summary', str(prefix)+'.tree.json', '--', sys.executable, BULGE/'measure.py',
                        '--metrics', str(prefix)+'.measurement.json', '--stdout', str(prefix)+'.paf',
                        '--stderr', str(prefix)+'.stderr', '--timeout', '300', '--', *call], stdout=log)
                output = Path(str(prefix)+'.paf')
                if set(parse_rows(output.read_text())) != set(rows0):
                    raise ValueError('scoring changed anchors')
                if mode == 'cfd':
                    if output.read_bytes() != base.read_bytes():
                        raise ValueError('CFD default output changed')
                else:
                    converted, counts, reasons, _, _ = audit(base, output)
                    annotations = sidecar_check(prefix/'scoring/annotations.tsv', old, guides, refs)
                    for row, item in zip(converted, annotations):
                        if row['status'] != item['crispr_bulge_reason'] or (row['cb'] is not None and row['cb'] != float(item['crispr_bulge'])):
                            raise ValueError('PAF/sidecar neural score mismatch')
                    tally = {'counts': counts, 'unavailable': reasons}
                    if statuses is not None and tally != statuses:
                        raise ValueError('mode/repetition availability changed')
                    statuses = tally
                    worker = json.loads((prefix/'scoring/metadata.json').read_text())['worker']
                    if worker and worker['model_load_count'] != 5:
                        raise ValueError('models loaded more than once')
                    if a.official_python:
                        for row, raw, item in zip(converted, old.splitlines(), annotations):
                            if row['cb'] is None:
                                continue
                            features, _ = independent_features(raw, guides, refs)
                            # One official reference evaluation per panel/eligible row, reused for both CLI modes/repetitions.
                            key = (label, row['id'])
                            if repetition == 0 and mode == 'crispr-bulge':
                                official_requests.append(dict(id=len(official_requests), **features))
                                comparisons.append((key, row['cb']))
                            else:
                                recorded = next(v for k, v in comparisons if k == key)
                                if abs(recorded-row['cb']) > 1e-6:
                                    raise ValueError('mode/repetition score changed')
                measure = json.loads(Path(str(prefix)+'.measurement.json').read_text())
                tree = json.loads(Path(str(prefix)+'.tree.json').read_text())
                result = dict(mode=mode, repetition=repetition, measurement=measure, tree=tree,
                              paf_sha256=digest(output), rows=len(rows0))
                results.append(result)
        panel = dict(label=label, bases=sum(map(len, refs.values())), contigs=len(refs), guides=len(guides),
            mismatches=mismatches, expected=len(expected) if expected is not None else len(parse_rows(retained)),
            expected_basis='independent direct-string forward oracle' if expected is not None else 'pinned previously independently validated 172-panel PAF',
            recovered=len(rows0), missing=0, unexpected=0, duplicates=0, cfd=cfd_counts,
            baseline_byte_identical=True, reference_sha256=digest(ref), guides_sha256=digest(guidefile),
            baseline_paf_sha256=digest(base), scoring=statuses, runs=results)
        summary['panels'].append(panel)
        if label == 'N':
            for mode in ('cfd', 'crispr-bulge', 'both'):
                runs = [r for r in results if r['mode'] == mode]
                summary['benchmark'].append(dict(panel=label, mode=mode, repetitions=len(runs),
                    median_wall_seconds=statistics.median(r['measurement']['wall_seconds'] for r in runs),
                    max_individual_rss_kib=max(r['measurement']['max_individual_child_rss_kib'] for r in runs),
                    max_sampled_tree_rss_kib=max(r['tree']['sampled_peak_rss_kib'] for r in runs)))
        print(label, panel['recovered'], statuses, flush=True)
    if a.official_python:
        request = a.out/'official-requests.jsonl'
        request.write_text(''.join(json.dumps(x)+'\n' for x in official_requests))
        with (a.out/'official.stdout').open('x') as stdout, (a.out/'official.stderr').open('x') as stderr:
            run([a.official_python, '-I', BULGE/'official_reference.py', '--upstream', a.official_root,
                '--input', request, '--output', a.out/'official-scores.jsonl', '--metrics', a.out/'official-metrics.json'],
                stdout=stdout, stderr=stderr)
        official = [json.loads(x) for x in (a.out/'official-scores.jsonl').read_text().splitlines()]
        if [x['id'] for x in official] != list(range(len(comparisons))):
            raise ValueError('official count/IDs changed')
        error = max(abs(score-item['score']) for (_, score), item in zip(comparisons, official))
        if error > 1e-6:
            raise ValueError('official prediction disagreement')
        summary['official'] = dict(rows=len(official), max_absolute_error=error, tolerance=1e-6,
            metrics=json.loads((a.out/'official-metrics.json').read_text()),
            requests_sha256=digest(request), predictions_sha256=digest(a.out/'official-scores.jsonl'))
    for path, sha in pinned.items():
        if digest(path) != sha:
            raise ValueError('tool/source changed during validation: '+path)
    summary['sha256'] = pinned
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
