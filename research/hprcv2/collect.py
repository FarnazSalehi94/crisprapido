#!/usr/bin/env python3
"""Collect checksum-committed pilots; never turn pending/missing results into zero."""
import argparse
import csv
import json
from pathlib import Path
import re
import zlib

from inventory import digest


def compressed_bytes(path):
    compressor = zlib.compressobj(6, zlib.DEFLATED, 31)
    count = 0
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            count += len(compressor.compress(block))
    return count+len(compressor.flush())


def only(paths, name):
    found = [Path(p) for p in paths if Path(p).name == name]
    if len(found) != 1:
        raise ValueError('expected one '+name)
    return found[0]


def load_task(root, project, task):
    stage = '01_search' if project == 'sbds' else '01_pilot'
    marker = root/stage/'HG00735_pat'/task/'SUCCESS.json'
    if not marker.is_file():
        return dict(status='PENDING_OR_INCOMPLETE', success_marker=str(marker))
    result = json.loads(marker.read_text())
    if result['status'] != 'COMPLETE' or result['monitor']['reason'] != 'completed':
        raise ValueError('invalid completion marker')
    files = result['outputs']
    for p, spec in files.items():
        if digest(p) != spec['sha256']:
            raise ValueError('output changed: '+p)
    qc = json.loads(only(files, 'qc-summary.json').read_text())
    search = json.loads(only(files, 'search-metrics.json').read_text())
    if search.get('reused_search'):
        signature = json.loads(only(files, 'SEARCH_COMPLETE.json').read_text())['search_signature']
        while search.get('reused_search'):
            original = Path(search['original_checkpoint'])
            if digest(original) != search['original_checkpoint_sha256']:
                raise ValueError('resumed checkpoint changed')
            checkpoint = json.loads(original.read_text())
            if checkpoint['search_signature'] != signature:
                raise ValueError('resumed search identity changed')
            search = checkpoint['metrics']
        log = (original.parent/'search.stderr').read_text()
    else:
        log = only(files, 'search.stderr').read_text()
    found = re.findall(r'^states=(\d+) search_wall_seconds=([0-9.eE+-]+)$', log, re.M)
    if len(found) != 1:
        raise ValueError('missing/duplicate native query timing')
    states, native_seconds = int(found[0][0]), float(found[0][1])
    scoring = json.loads(only(files, 'scoring-metrics.json').read_text())
    ann = only(files, 'annotations.tsv')
    model_meta = json.loads((ann.parent/'metadata.json').read_text())
    if qc['cb_statuses'].get('scored', 0) and model_meta['worker']['model_load_count'] != 5:
        raise ValueError('ensemble not loaded exactly once')
    enriched = dict(status='COMPLETE', completion_marker=str(marker), completion_sha256=digest(marker),
                    qc=qc, search_metrics=search, scoring_metrics=scoring,
                    native_states=states, native_search_seconds_including_replay=native_seconds,
                    native_batch_mean_seconds_per_guide=native_seconds/len(qc['per_guide']),
                    other_search_wall_residual_seconds=search['wall_seconds']-native_seconds,
                    residual_scope='Everything else in search call: validation/staging, loads, Python/Rust, CFD and output; not isolated index time.',
                    monitor=result['monitor'], scorer_worker=model_meta['worker'],
                    retained_artifact_bytes=sum(s['bytes'] for s in files.values()),
                    full_paf_gzip_bytes=compressed_bytes(only(files, 'eligible.both.paf')),
                    annotations_gzip_bytes=compressed_bytes(ann),
                    original_cfd_paf_gzip_bytes=compressed_bytes(only(files, 'eligible.cfd.paf')))
    return enriched


def site_records(path):
    result = {}
    with Path(path).open() as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            tags = {s.split(':', 2)[0]: s.split(':', 2)[2] for s in fields[12:]}
            key = (fields[0], fields[5], fields[4], fields[8] if fields[4] == '+' else fields[7])
            if key in result:
                raise ValueError('duplicate site')
            result[key] = dict(span=fields[7:9], cigar=tags['cg'], cfd=tags.get('cf'), cb=tags.get('cb'), status=tags.get('bu'))
    return result


def precision_views(files):
    scores = {}
    with only(files, 'records.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = (row['guide'], row['contig'], row['strand'], row['end0'] if row['strand'] == '+' else row['start0'])
            scores[key] = None if row['independent_cfd'] in ('', 'None') else float(row['independent_cfd'])
    return {str(threshold): {k for k, value in scores.items() if value is None or value >= threshold}
            for threshold in [.05, .20]}


def concordance(first, second):
    x = json.loads(Path(first['completion_marker']).read_text())
    y = json.loads(Path(second['completion_marker']).read_text())
    a = site_records(only(x['outputs'], 'eligible.both.paf'))
    b = site_records(only(y['outputs'], 'eligible.both.paf'))
    shared = a.keys() & b.keys()
    av, bv = precision_views(x['outputs']), precision_views(y['outputs'])
    return dict(full_precision_cfd_view_concordance={threshold: dict(first_retained=len(av[threshold]), second_retained=len(bv[threshold]),
                    first_only=len(av[threshold]-bv[threshold]), second_only=len(bv[threshold]-av[threshold])) for threshold in av},
                shared_neural_available_to_unavailable=sum(a[k]['cb'] is not None and b[k]['cb'] is None for k in shared),
                first_sites=len(a), second_sites=len(b), shared=len(shared), first_only=len(a.keys()-b.keys()),
                second_only=len(b.keys()-a.keys()),
                shared_representative_differences=sum((a[k]['span'], a[k]['cigar']) != (b[k]['span'], b[k]['cigar']) for k in shared),
                shared_stored_cfd_differences=sum(a[k]['cfd'] != b[k]['cfd'] for k in shared),
                caveat='Exact same-guide/contig/strand/PAM-boundary keys. Stored CFD is rounded. Changes in representative can change score; no biological preference inferred.')


def collect(root, project):
    names = ['primary', 'historical-domain-control'] if project == 'sbds' else ['A', 'B']
    results = {name: load_task(root, project, name) for name in names}
    complete = all(v['status'] == 'COMPLETE' for v in results.values())
    report = dict(project=project, status='COMPLETE_SAFETY_PILOT' if complete else 'PENDING_OR_INCOMPLETE',
                  scope='One assembly, not a full-cohort result or validated library-wide cost model.',
                  full_brunello_authorized=False, configurations=results)
    if complete:
        first, second = ('historical-domain-control', 'primary') if project == 'sbds' else ('A', 'B')
        report['nested_domain_concordance'] = concordance(results[first], results[second])
        if report['nested_domain_concordance']['first_only']:
            report['expansion_gate'] = 'BLOCKED: narrower-domain sites missing in broader domain; investigate'
        else:
            report['expansion_gate'] = 'Still blocked pending input integrity, broader sampling and explicit full-Brunello review'
        if project == 'brunello':
            combinations = 77441*464
            report['unvalidated_linear_sensitivity_scenarios'] = {}
            for name, value in results.items():
                guide_count = len(value['qc']['per_guide'])
                report['unvalidated_linear_sensitivity_scenarios'][name] = dict(
                    guide_haplotype_combinations=combinations,
                    native_single_core_hours=value['native_batch_mean_seconds_per_guide']*combinations/3600,
                    full_eligible_paf_gzip_bytes=value['full_paf_gzip_bytes']/guide_count*combinations,
                    full_annotations_gzip_bytes=value['annotations_gzip_bytes']/guide_count*combinations,
                    caveat='Three guides, one assembly, one fresh call/config, concurrent pilot neighbour and uncontrolled caches. NOT a confidence interval, resource reservation or launch approval. Excludes batching/setup, scoring, validation, failures and mapping costs. No individual per-guide runtime instrumentation in this attested binary.')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', required=True, choices=['sbds', 'brunello'])
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    if a.root.name != a.project+'_pam_aware_v1':
        raise ValueError('cross-project root')
    report = collect(a.root, a.project)
    with a.out.open('x') as f:
        json.dump(report, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(status=report['status'], output=str(a.out))))


if __name__ == '__main__':
    main()
