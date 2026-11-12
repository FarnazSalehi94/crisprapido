#!/usr/bin/env python3
"""Analyze the completed 16-guide benchmark: configuration comparison and nesting.

Emitted-witness comparisons are bounded to this panel/assemblies; nothing here
is a genome-wide recall claim. Full-precision CFD comes from the independent
audit, never rounded PAF tags.
"""
import argparse
import collections
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
for candidate in (HERE, HERE.parent, HERE.parent.parent):
    sys.path.insert(0, str(candidate))
from collect import compressed_bytes, only  # noqa: E402
from inventory import digest  # noqa: E402

CHAIN = ['historical', 'intermediate', 'candidate-A', 'candidate-B']
CFD_BINS = [(.02, 'ge_0.02'), (.05, 'ge_0.05'), (.10, 'ge_0.10'),
            (.20, 'ge_0.20'), (.30, 'ge_0.30'), (.50, 'ge_0.50')]


def anchor_key(fields):
    return (fields[0], fields[5], fields[4], fields[8] if fields[4] == '+' else fields[7])


def site_records(path):
    """Anchor key -> (span, cigar, cfd, cb, status)."""
    result = {}
    with Path(path).open() as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            tags = {s.split(':', 2)[0]: s.split(':', 2)[2] for s in fields[12:]}
            key = anchor_key(fields)
            if key in result:
                raise ValueError('duplicate anchor: %r' % (key,))
            result[key] = dict(span=fields[7:9], cigar=tags['cg'], cfd=tags.get('cf'),
                               cb=tags.get('cb'), status=tags.get('bu'))
    return result


def audited_cfd(records_path):
    """Anchor key -> full-precision independently reconstructed CFD (or None)."""
    import csv
    result = {}
    with Path(records_path).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = (row['guide'], row['contig'], row['strand'],
                   row['end0'] if row['strand'] == '+' else row['start0'])
            value = row['independent_cfd']
            result[key] = None if value in ('', 'None') else float(value)
    return result


def cfd_profile(values):
    defined = [v for v in values if v is not None]
    profile = dict(total=len(values), defined=len(defined), unavailable=len(values)-len(defined))
    for threshold, name in CFD_BINS:
        profile[name] = sum(v >= threshold for v in defined)
    if defined:
        ordered = sorted(defined)
        profile['median'] = ordered[len(ordered)//2]
        profile['max'] = ordered[-1]
    return profile


def nested_comparison(sets):
    """Per-step additions and shared-representative changes along the config chain."""
    result = {}
    for narrower, broader in zip(CHAIN, CHAIN[1:]):
        if narrower not in sets or broader not in sets:
            continue
        small, large = sets[narrower], sets[broader]
        shared = set(small) & set(large)
        added = set(large) - set(small)
        result[(narrower, broader)] = dict(
            narrower_sites=len(small), broader_sites=len(large), shared=len(shared),
            added_sites=len(added), removed_sites=len(set(small)-set(large)),
            shared_representative_differences=sum(
                (small[k]['span'], small[k]['cigar']) != (large[k]['span'], large[k]['cigar']) for k in shared),
            shared_stored_cfd_differences=sum(
                small[k]['cfd'] is not None and large[k]['cfd'] is not None
                and float(small[k]['cfd']) != float(large[k]['cfd']) for k in shared),
            shared_neural_availability_changes=sum(
                (small[k]['cb'] is not None) != (large[k]['cb'] is not None) for k in shared),
            caveat='Exact same-guide/contig/strand/PAM-boundary keys; representative policy is '
                   'edit/mismatch/CIGAR/start-ordered, not CFD-ordered; no biological preference.')
    return result


def load_run(root, assembly, config):
    marker = Path(root)/'01_pilot'/'expanded-benchmark-v1'/assembly/config/'SUCCESS.json'
    if not marker.is_file():
        return dict(status='PENDING_OR_INCOMPLETE')
    result = json.loads(marker.read_text())
    for path, spec in result['outputs'].items():
        require_ok = digest(path) == spec['sha256']
        if not require_ok:
            raise ValueError('completed output changed: '+path)
    qc = json.loads(only(result['outputs'], 'qc-summary.json').read_text())
    search = json.loads(only(result['outputs'], 'search-metrics.json').read_text())
    scoring = json.loads(only(result['outputs'], 'scoring-metrics.json').read_text())
    audit = json.loads(only(result['outputs'], 'audit-metrics.json').read_text())
    paf = only(result['outputs'], 'eligible.both.paf')
    annotations = only(result['outputs'], 'annotations.tsv')
    records = only(result['outputs'], 'records.tsv')
    return dict(status='COMPLETE', qc=qc, search=search, scoring=scoring, audit=audit,
                monitor=result['monitor'], completion_sha256=digest(marker),
                sites=site_records(paf), audited_cfd=audited_cfd(records),
                full_paf_gzip_bytes=compressed_bytes(paf),
                annotations_gzip_bytes=compressed_bytes(annotations),
                raw_paf_bytes=Path(paf).stat().st_size,
                raw_annotations_bytes=Path(annotations).stat().st_size)


def collect(root, out):
    runs = {assembly: {config: load_run(root, assembly, config)
                       for config in CHAIN}
            for assembly in ['HG00735_pat', 'HG01496_pat']}
    if any(run['status'] != 'COMPLETE' for a in runs for run in runs[a].values()):
        return dict(status='PENDING_OR_INCOMPLETE', runs={a: {c: r['status'] for c, r in v.items()} for a, v in runs.items()})
    report = dict(status='COMPLETE', scope='16-guide panel, two assemblies, four configurations')
    per_run = {}
    for assembly, by_config in runs.items():
        for config, run in by_config.items():
            qc = run['qc']
            per_guide_sites = {name: stats['eligible_hits'] for name, stats in qc['per_guide'].items()}
            additional = None
            if config != CHAIN[0]:
                narrower = CHAIN[CHAIN.index(config)-1]
                shared = set(runs[assembly][narrower]['sites']) & set(run['sites'])
                added = set(run['sites']) - set(runs[assembly][narrower]['sites'])
                additional = dict(count=len(added), cfd=cfd_profile([run['audited_cfd'][k] for k in added]),
                                  neural_available=sum(run['sites'][k]['cb'] is not None for k in added))
            per_run[(assembly, config)] = dict(
                search_wall_seconds=run['search']['wall_seconds'],
                native_traversal_seconds=run['qc']['native_traversal_seconds'],
                native_mean_seconds_per_guide=run['qc']['native_traversal_seconds']/16,
                native_states=run['qc']['native_states'],
                audit_wall_seconds=run['audit']['wall_seconds'],
                scoring_wall_seconds=run['scoring']['wall_seconds'],
                workflow_wall_seconds=run['monitor']['elapsed_seconds'],
                peak_tree_rss_kib=run['monitor']['sampled_peak_rss_kib'],
                rows=qc['eligible_hits'], cb_statuses=qc['cb_statuses'],
                per_guide_sites=per_guide_sites,
                per_guide_site_spread=dict(min=min(per_guide_sites.values()), max=max(per_guide_sites.values())),
                all_sites_cfd=cfd_profile(list(run['audited_cfd'].values())),
                views={v: sum(s['view_'+v+'_rows'] for s in qc['per_guide'].values())
                       for v in ['none', '0.05', '0.20']},
                view_gzip_bytes=qc['analytical_view_gzip_bytes'],
                full_paf_gzip_bytes=run['full_paf_gzip_bytes'],
                annotations_gzip_bytes=run['annotations_gzip_bytes'],
                raw_paf_bytes=run['raw_paf_bytes'], raw_annotations_bytes=run['raw_annotations_bytes'],
                additional_vs_narrower=additional)
    report['runs'] = {f'{a}/{c}': v for (a, c), v in per_run.items()}
    report['nested'] = {assembly: {f'{n}->{b}': v for (n, b), v in nested_comparison(
        {c: run['sites'] for c, run in runs[assembly].items()}).items()} for assembly in runs}
    report['caveats'] = [
        'Single fresh call per (assembly, config); concurrent sibling chain on the same node; uncontrolled caches.',
        '16 timing-blind, deterministically selected guides; not a random library sample.',
        'Additional-site CFD profiles use independent full-precision reconstruction; unavailable CFD is kept separate.',
        'Nested set comparisons are exact-anchor keys; representative changes can alter scores without search loss.',
        'No genome-wide recall, biological calibration or full-library timing claim.']
    Path(out).write_text(json.dumps(report, indent=2, default=str)+'\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    result = collect(a.root, a.out)
    print(json.dumps(dict(status=result['status'], output=str(a.out))))


if __name__ == '__main__':
    main()
