#!/usr/bin/env python3
"""Deterministic 16-guide panel from historical retained counts + GC strata.

Historical CFD>=0.20 retained counts are a biased off-target proxy (thresholded,
legacy PAM semantics); the panel supports workload/biology comparison across
configurations, never a population estimate. Guides are not cherry-picked by
measured timing: within each stratum cell the median-count guide is taken with
an identifier tie-break.
"""
import argparse
import collections
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from inventory import LIBRARY_SHA256, digest, library


def historical_counts(zst, assembly):
    proc = subprocess.Popen(['zstd', '-dc', str(zst)], stdout=subprocess.PIPE, text=True)
    counts = collections.Counter()
    for line in proc.stdout:
        counts[line.split('\t', 1)[0]] += 1
    if proc.wait() != 0 or not counts:
        raise ValueError('historical retained stream unreadable: '+assembly)
    return counts


def bounds_of(values):
    nonzero = sorted(v for v in values if v)
    if not nonzero:
        raise ValueError('no nonzero counts')
    def pct(p):
        return nonzero[min(len(nonzero)-1, int(p*len(nonzero)))]
    return dict(zero=(0, 1), low=(1, pct(.25)+1), medium=(pct(.25)+1, pct(.75)+1),
                high=(pct(.75)+1, pct(.90)+1), top=(pct(.90)+1, max(nonzero)+1))


def select(records, counts, assemblies):
    gc = {name: sum(c in 'GC' for c in seq) for name, seq in records}
    gc_values = sorted(gc.values())
    third, two_thirds = gc_values[len(gc_values)//3], gc_values[2*len(gc_values)//3]
    def gc_stratum(value):
        return 'low' if value <= third else 'high' if value >= two_thirds else 'mid'
    total = {name: sum(counts[a][name] for a in assemblies) for name, _ in records}
    bounds = bounds_of(total.values())
    def count_stratum(value):
        for name in ['zero', 'low', 'medium', 'high', 'top']:
            if bounds[name][0] <= value < bounds[name][1]:
                return name
        raise ValueError('unstratified count: %d' % value)
    cells = collections.defaultdict(list)
    for name, sequence in records:
        if 'Non-Targeting_Control' in name:
            continue
        cells[(count_stratum(total[name]), gc_stratum(gc[name]))].append((total[name], name, sequence))
    panel, used = [], set()
    def take(group, offset):
        ordered = sorted(group)
        index = min(len(ordered)-1, len(ordered)//2+offset)
        for _, name, sequence in ordered[index:]:
            if name not in used:
                return name, sequence
        return None, None
    for key in sorted(cells):
        name, sequence = take(cells[key], 0)
        if name is None:
            continue
        used.add(name)
        panel.append(dict(guide_id=name, sequence=sequence,
                          stratum='count:'+key[0]+'_gc:'+key[1], gc_bases=gc[name],
                          historical_retained={a: counts[a][name] for a in assemblies},
                          historical_retained_total=total[name]))
    controls = [name for name, _ in records if 'Non-Targeting_Control' in name]
    if controls:
        chosen = min(controls, key=lambda n: (total[n], n))
        sequence = dict(records)[chosen]
        if chosen not in used and len(panel) < 16:
            used.add(chosen)
            panel.append(dict(guide_id=chosen, sequence=sequence, stratum='non_targeting_control',
                              gc_bases=gc[chosen],
                              historical_retained={a: counts[a][chosen] for a in assemblies},
                              historical_retained_total=total[chosen]))
    while len(panel) < 16:
        progressed = False
        for key in sorted(cells, key=lambda k: -len(cells[k])):
            if len(panel) >= 16:
                break
            for offset in range(1, 8):
                name, sequence = take(cells[key], len(used & {n for _, n, _ in cells[key]}))
                if name is None or name in used:
                    continue
                used.add(name)
                panel.append(dict(guide_id=name, sequence=sequence,
                                  stratum='count:'+key[0]+'_gc:'+key[1]+'_fill', gc_bases=gc[name],
                                  historical_retained={a: counts[a][name] for a in assemblies},
                                  historical_retained_total=total[name]))
                progressed = True
                break
            if len(panel) >= 16:
                break
        if not progressed:
            raise ValueError('cannot fill panel to 16 from strata')
    panel = panel[:16]
    if len({r['guide_id'] for r in panel}) != len(panel) or len(panel) < 16:
        raise ValueError('panel uniqueness/size failure: %d' % len(panel))
    return panel


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library', required=True, type=Path)
    p.add_argument('--assembly-one', required=True)
    p.add_argument('--assembly-two', required=True)
    p.add_argument('--retained-one', required=True, type=Path)
    p.add_argument('--retained-two', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    if digest(a.library) != LIBRARY_SHA256:
        raise ValueError('frozen library checksum mismatch')
    assemblies = [a.assembly_one, a.assembly_two]
    records = library(a.library)
    counts = {assemblies[0]: historical_counts(a.retained_one, a.assembly_one),
              assemblies[1]: historical_counts(a.retained_two, a.assembly_two)}
    panel = select(records, counts, assemblies)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('x') as f:
        for r in panel:
            f.write('\t'.join([r['guide_id'], r['sequence'], r['stratum'], str(r['gc_bases']),
                               json.dumps(r['historical_retained']), str(r['historical_retained_total'])])+'\n')
    fasta = a.out.with_suffix('.fa')
    with fasta.open('x') as f:
        for r in panel:
            f.write('>'+r['guide_id']+'\n'+r['sequence']+'\n')
    print(json.dumps(dict(panel=str(a.out), fasta=str(fasta), sha256=digest(fasta),
                          guides=[r['guide_id'] for r in panel],
                          strata=sorted({r['stratum'] for r in panel}),
                          caveat='Median-count guide per count/GC cell; historical retained-count proxy is thresholded and from legacy-semantics runs; includes at most one non-targeting control.'), indent=2))


if __name__ == '__main__':
    main()
