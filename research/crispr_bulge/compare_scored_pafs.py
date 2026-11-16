#!/usr/bin/env python3
"""Validate chunked-scoring equivalence against the single-pass reference.

Chunked scoring changes model-worker batch composition (one partial batch per
chunk instead of one per file). TensorFlow CPU inference is batch-composition
sensitive in low-order ULPs, so a small fraction of near-zero scores can
differ in trailing digits. This tool makes the acceptance rule explicit:
rows must be identical in order and content except the cb:f tag, which may
differ by at most TOLERANCE absolute (default 1e-9; observed worst case
~1e-14). Any other difference fails validation.
"""
import argparse
import json
import math
from pathlib import Path

TOLERANCE = 1e-9


def tags(line):
    fields = line.rstrip('\n').split('\t')
    return fields[:12], {f.split(':', 2)[0]: f for f in fields[12:]}


def compare(chunked, reference):
    differences = []
    max_cb_delta = 0.0
    rows = 0
    with open(chunked) as a, open(reference) as b:
        for n, (x, y) in enumerate(zip(a, b), 1):
            rows += 1
            if x == y:
                continue
            fx, tx = tags(x)
            fy, ty = tags(y)
            if fx != fy or set(tx) != set(ty):
                differences.append((n, 'structure', x, y))
                continue
            for key in tx:
                if key == 'cb':
                    vx = float(tx[key].rsplit(':', 1)[1])
                    vy = float(ty[key].rsplit(':', 1)[1])
                    delta = abs(vx-vy)
                    max_cb_delta = max(max_cb_delta, delta)
                    if delta > TOLERANCE or not (math.isfinite(vx) and math.isfinite(vy)):
                        differences.append((n, 'cb '+key, x, y))
                elif tx[key] != ty[key]:
                    differences.append((n, key, x, y))
    if a_extra := sum(1 for _ in open(chunked)) - rows:
        differences.append((-1, 'chunked row count', '', ''))
    if b_extra := sum(1 for _ in open(reference)) - rows:
        differences.append((-2, 'reference row count', '', ''))
    return rows, differences, max_cb_delta


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--chunked', required=True)
    p.add_argument('--reference', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    rows, differences, max_cb_delta = compare(a.chunked, a.reference)
    report = dict(rows=rows, differences=len(differences),
                  max_cb_absolute_delta=max_cb_delta, tolerance=TOLERANCE,
                  chunked=a.chunked, reference=a.reference,
                  status='PASS' if not differences else 'FAIL',
                  note='identical rows/order; cb:f differences bounded by batch-'
                       'composition ULP effects at near-zero scores')
    Path(a.output).write_text(json.dumps(report, indent=2)+'\n')
    if differences:
        for n, kind, x, y in differences[:5]:
            print('row %s: %s\n chunked: %s\n reference: %s' % (n, kind, x.strip(), y.strip()))
    print(json.dumps({k: report[k] for k in ('status', 'rows', 'differences', 'max_cb_absolute_delta')}))
    return 0 if not differences else 1


if __name__ == '__main__':
    raise SystemExit(main())
