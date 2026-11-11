#!/usr/bin/env python3
"""Aggregate the bounded benchmark; explicit unvalidated linear extrapolations only."""
import argparse
import json
from pathlib import Path

from bench_worker import parse_time_v


def extrapolate(points, total):
    """Linear projection through the two largest measured points; no uncertainty claim."""
    ordered = sorted(points)
    (x1, y1), (x2, y2) = ordered[-2], ordered[-1]
    if x2 == x1:
        return None
    return y1 + (y2-y1)*(total-x1)/(x2-x1)


def load_runs(root, engine, sizes):
    stage = Path(root)/('01_new_backend' if engine == 'new' else '02_old_pipeline')
    runs = {}
    for size in sizes:
        marker = stage/('N%d' % size)/'SUCCESS.json'
        if not marker.is_file():
            runs[size] = dict(status='PENDING_OR_INCOMPLETE')
            continue
        result = json.loads(marker.read_text())
        runs[size] = dict(status='COMPLETE', summary=json.loads(next(stage.glob('N%d/attempt-*/summary.json' % size)).read_text()),
                          monitor=result['monitor'])
    return runs


def historical(base):
    base = Path(base)
    brunello = parse_time_v((base/'results/Brunello_all_HPRCv2_scored_unscored_v2/HG00735_pat/brunello.time').read_text())
    meta = dict(line.split('\t', 1) for line in (base/'results/Brunello_all_HPRCv2_scored_unscored_v2/HG00735_pat/brunello.metadata.tsv').read_text().splitlines())
    brunello.update(guides=int(meta['guide_count']), allocated_cpus=int(meta['allocated_cpus']),
                    application_threads=int(meta['application_threads']), node=meta['node'],
                    total_candidates=int(meta['total_candidates']),
                    description='Preserved production full-library run, one assembly, m2/b1/z1/f0.75/GG, index pre-staged on local scratch')
    sbds = parse_time_v((base/'results/SBDS_all_HPRCv2/HG00735_pat/SBDS.time').read_text())
    sbds.update(guides=1, description='Preserved production single-guide SBDS run, same assembly, m2/b1/z1, index read from MooseFS, -t 1')
    return dict(full_library=brunello, sbds_single=sbds)


def collect(root, base, out):
    sizes = [1, 100, 1000]
    engines = {name: load_runs(root, name, sizes) for name in ['new', 'old']}
    if any(r['status'] != 'COMPLETE' for e in engines.values() for r in e.values()):
        return dict(status='PENDING_OR_INCOMPLETE', engines=engines)
    total_guides = 77441
    result = dict(status='COMPLETE_BOUNDED_BENCHMARK', total_guides=total_guides, sizes=sizes)
    walls = {name: {size: r['summary']['search_wall_seconds'] if name == 'new' else r['summary']['wall_seconds']
                    for size, r in runs.items()} for name, runs in engines.items()}
    for name in ['new', 'old']:
        runs = engines[name]
        points = [(size, walls[name][size]) for size in sizes]
        est = extrapolate(points, total_guides)
        entry = dict(runs={size: runs[size]['summary'] for size in sizes},
                     measured_wall_seconds=walls[name],
                     fixed_overhead_single_guide_seconds=walls[name][1]-engines[name][1]['summary'].get('native_traversal_seconds', 0) if name == 'new' else walls[name][1],
                     marginal_wall_seconds_per_guide=(walls[name][1000]-walls[name][100])/900,
                     linear_extrapolated_full_library_seconds=est,
                     linear_extrapolated_full_library_hours=est/3600,
                     extrapolated_throughput_guides_per_hour=total_guides/(est/3600) if est else None)
        if name == 'new':
            tr = {size: runs[size]['summary']['native_traversal_seconds'] for size in sizes}
            rows = {size: runs[size]['summary']['rows_emitted'] for size in sizes}
            audit = {size: runs[size]['summary']['audit_wall_seconds'] for size in sizes}
            entry.update(native_traversal_seconds=tr,
                         marginal_traversal_seconds_per_guide=(tr[1000]-tr[100])/900,
                         rows_emitted=rows,
                         extrapolated_rows_full_library=rows[1000]/1000*total_guides,
                         extrapolated_audit_seconds_full_library=audit[1000]/rows[1000]*(rows[1000]/1000*total_guides),
                         native_states={size: runs[size]['summary']['native_states'] for size in sizes})
        else:
            entry.update(historical_validation=dict(
                extrapolated_seconds=est,
                actual_measured_full_library_seconds=None,  # filled by caller against preserved run
                note='Compare with historical full-library measurement below; closeness validates the linear model, not a guarantee.'))
        result[name] = entry
    hist = historical(base)
    result['old']['historical_validation']['actual_measured_full_library_seconds'] = hist['full_library']['wall_seconds']
    result['old']['historical_validation']['extrapolation_vs_actual_ratio'] = (
        result['old']['linear_extrapolated_full_library_seconds']/hist['full_library']['wall_seconds'])
    result['historical'] = hist
    result['historical']['full_library']['throughput_guides_per_hour'] = total_guides/(hist['full_library']['wall_seconds']/3600)
    result['cpu_utilization'] = dict(
        new='1 of 8 allocated CPUs used; experimental backend is single-threaded by design and rejects -t != 1',
        old='measured cpu_percent per run; historical full library averaged 190% of one CPU despite -t 16')
    result['full_hprcv2_464_assemblies_sequential_hours'] = {
        name: result[name]['linear_extrapolated_full_library_hours']*464 for name in ['new', 'old']}
    result['caveats'] = [
        'First 1/100/1000 guides of the frozen library on ONE verified assembly; not a random or complete sample.',
        'Single fresh call per (engine, size) on a quiet sequential node; uncontrolled caches; not repeated medians.',
        'Linear extrapolation through the two largest measured points only; validated only where the preserved full-library measurement exists.',
        'The new engine runs strict NGG with corrected eligibility; the old engine emits legacy non-NGG/geometry rows, so row counts are NOT directly comparable.',
        'Per-assembly one-time costs (index build/import/verification, staging) are excluded from marginal rates and listed separately.',
        'No full-library or cohort launch is authorized by these numbers.']
    Path(out).write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--base', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    result = collect(a.root, a.base, a.out)
    print(json.dumps(dict(status=result['status'], output=str(a.out))))


if __name__ == '__main__':
    main()
