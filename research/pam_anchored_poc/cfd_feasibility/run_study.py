#!/usr/bin/env python3
"""Bounded offline comparisons; runs unchanged native backend only on tiny fixtures."""
import argparse
import csv
import hashlib
import json
import platform
import resource
import subprocess
import tempfile
import time
from pathlib import Path

from study import (ROOT, fixtures, keep, load_tables, numeric_upper, oracle,
                   representatives, reverse_model, threshold)


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def native_check(binary, index, refpath, guidepath, config, expected):
    m, z, frac, pam = config
    run = subprocess.run([binary, str(index), str(refpath), str(guidepath), str(m), str(z),
                          str(frac), pam, str(int(z > 0))], text=True, capture_output=True,
                         check=True, timeout=30)
    lines = run.stdout.splitlines()
    start = next(i for i, s in enumerate(lines) if s.startswith('guide\tcontig\t'))
    got = {}
    for row in csv.DictReader(lines[start:], delimiter='\t'):
        key = row['contig'], row['strand'], int(row['anchor'])
        assert key not in got
        got[key] = int(row['start']), int(row['end']), row['cigar'], int(row['mismatches'])
    want = {k: (a.start, a.end, a.cigar, a.mm) for k, a in expected.items()}
    assert got == want, (config, len(got), len(want), list(set(want)-set(got))[:3])
    return dict(records=len(got), matching_representatives=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True, type=Path, help='new exclusive small-results directory')
    p.add_argument('--builder', required=True)
    p.add_argument('--native', required=True)
    p.add_argument('--cfd-probe', required=True)
    a = p.parse_args()
    # Risk reduction for this small local study; not a process-tree memory cap.
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    a.out.mkdir(exist_ok=False)
    tables = load_tables()
    guide, resolved, nrefs = fixtures()
    configurations = [(5, 0, .75, 'GG'), (5, 1, .75, 'GG'), (5, 2, .75, 'GG'),
                      (8, 2, .5, 'GG'), (20, 0, 0, 'GG')]
    configs = configurations + [(5, 0, .75, x+y) for x in 'ACGT' for y in 'ACGT' if x+y != 'GG']
    result = dict(scope='synthetic bounded occurrence-set model, not pruned FM-index performance',
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  python=platform.python_version(), host=platform.node(), guide=guide,
                  resolved_contigs=len(resolved), n_panel_contigs=len(nrefs),
                  resolved_bases=sum(map(len, resolved.values())), n_panel_bases=sum(map(len, nrefs.values())),
                  configurations=[], representative_counterexamples=[])
    result['sha256'] = {str(path): digest(path) for path in (
        ROOT/'mismatch_scores.txt', ROOT/'pam_scores.txt', ROOT/'src/cfd_score.rs',
        ROOT/'research/pam_anchored_poc/poc.cpp', ROOT/'scripts/release_validation/audit_emitted.py',
        Path(__file__), Path(__file__).with_name('study.py'),
        Path(a.native), Path(a.builder), Path(a.cfd_probe))}
    from sys import path
    path.insert(0, str(ROOT/'research/pam_anchored_poc/compatibility'))
    from index_manifest import bind_build
    probes = {}
    start_total = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='crispr-cfd-feasibility-') as tmp:
        tmp = Path(tmp)
        guidepath = tmp/'guide.fa'
        guidepath.write_text('>g\n'+guide+'\n')
        for label, refs in [('resolved', resolved), ('N', nrefs)]:
            refpath, index = tmp/(label+'.fa'), tmp/(label+'_index')
            refpath.write_text(''.join(f'>{c}\n{s}\n' for c, s in refs.items()))
            # Only tiny synthetic indexes are built; HPRC indexes are untouched.
            with (a.out/(label+'_build.log')).open('w') as log:
                subprocess.run(['python3', str(ROOT/'research/pam_anchored_poc/compatibility/index_manifest.py'),
                                'build', '--builder', a.builder, '--reference', str(refpath), '--index', str(index)],
                               stdout=log, stderr=subprocess.STDOUT, check=True, timeout=30)
            result[label+'_fasta_sha256'] = digest(refpath)
            for config in configs if label == 'resolved' else configurations:
                m, z, frac, pam = config
                started = time.perf_counter()
                truth = oracle(refs, guide, pam, m, z, frac, tables)
                oracle_seconds = time.perf_counter()-started
                canonical = representatives(truth)
                unpruned, base_stats = reverse_model(refs, guide, pam, m, z, frac, tables)
                assert unpruned == truth, (label, config, len(unpruned), len(truth))
                native = native_check(a.native, index, refpath, guidepath, config, canonical)
                for hit in truth:
                    probes[(guide, hit.target, hit.cigar, pam)] = hit.score
                row = dict(panel=label, m=m, z=z, fraction=frac, pam=pam,
                           oracle_alignment_paths=len(truth), oracle_anchors=len(canonical),
                           unavailable_representatives=sum(x.score is None for x in canonical.values()),
                           oracle_seconds=oracle_seconds, unpruned_stats=base_stats, native=native, thresholds=[])
                for tau in (0.0, .01, .1, .5, 1.0):
                    timings = []
                    # Alternating repeated model timings; refinement cost below is separate.
                    for _ in range(3):
                        t = time.perf_counter()
                        base, bs = reverse_model(refs, guide, pam, m, z, frac, tables)
                        base_time = time.perf_counter()-t
                        t = time.perf_counter()
                        discovered, stats = reverse_model(refs, guide, pam, m, z, frac, tables, tau)
                        timings.append(dict(unpruned_seconds=base_time, discovery_seconds=time.perf_counter()-t))
                    wanted_paths = {x for x in truth if keep(x, tau)}
                    assert wanted_paths <= discovered, (label, config, tau, 'lost supported/unknown path')
                    keys = {x.key for x in discovered}
                    # Research-only refinement using independently enumerated full paths.
                    # Not a measured implementation of a future targeted anchor verifier.
                    t = time.perf_counter()
                    refined = representatives(x for x in truth if x.key in keys)
                    final = {k: v for k, v in refined.items() if keep(v, tau)}
                    refine_lookup_seconds = time.perf_counter()-t
                    want = {k: v for k, v in canonical.items() if keep(v, tau)}
                    assert final == want
                    naive = {k: v for k, v in representatives(discovered).items() if keep(v, tau)}
                    diffs = [k for k in set(naive) | set(want) if naive.get(k) != want.get(k)]
                    if diffs and len(result['representative_counterexamples']) < 3:
                        k = diffs[0]
                        def describe(hit):
                            return None if hit is None else dict(cigar=hit.cigar, score=hit.score,
                                                                start=hit.start, end=hit.end, order=hit.order())
                        result['representative_counterexamples'].append(dict(panel=label, config=config, threshold=tau,
                            key=k, canonical_before_filter=describe(canonical.get(k)),
                            expected_after_filter=describe(want.get(k)), naive_after_pruning=describe(naive.get(k))))
                    row['thresholds'].append(dict(threshold=tau, final_anchors=len(final),
                        retained_unavailable=sum(x.score is None for x in final.values()),
                        discovery_anchors=len(keys), stats=stats,
                        state_reduction_fraction=1-stats['states']/bs['states'],
                        missing_paths=0, missing_or_unexpected_final_anchors=0, representative_differences_after_refinement=0,
                        naive_representative_or_membership_differences=len(diffs), timings=timings,
                        refinement_from_precomputed_oracle_seconds=refine_lookup_seconds))
                result['configurations'].append(row)
                print(label, config, 'anchors', len(canonical), 'states', base_stats['states'], flush=True)
        # Compare every distinct oracle reconstruction against the unchanged Rust scorer.
        inputs = sorted(probes)
        run = subprocess.run([a.cfd_probe], input=''.join('\t'.join(x)+'\n' for x in inputs), cwd=ROOT,
                             capture_output=True, text=True, check=True, timeout=30)
        actual = run.stdout.splitlines()
        assert len(actual) == len(inputs)
        for fields, value in zip(inputs, actual):
            expected = probes[fields]
            assert ((value == 'NA' and expected is None) or
                    (value != 'NA' and expected is not None and abs(float(value)-expected) < 1e-9)), (fields, value, expected)
        result['rust_cfd_comparison'] = dict(unique_inputs=len(inputs),
            defined=sum(x is not None for x in probes.values()), unavailable=sum(x is None for x in probes.values()), differences=0)
    result['elapsed_seconds'] = time.perf_counter()-start_total
    (a.out/'summary.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('PASS', result['rust_cfd_comparison'], 'elapsed', result['elapsed_seconds'])


if __name__ == '__main__':
    main()
