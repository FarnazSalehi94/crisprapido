#!/usr/bin/env python3
"""Verify full-run completeness across the assembly manifest; never zero-fill.

A missing SUCCESS marker, a changed checksum, a failed audit or an unparsable
summary leaves the assembly marked INCOMPLETE in the ledger. Absence is never
reported as a completed empty search.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
for candidate in (HERE, HERE.parent, HERE.parent.parent):
    sys.path.insert(0, str(candidate))
from inventory import digest  # noqa: E402


def verify(root, manifest, plan, out):
    assemblies = [line.split('\t', 1) for line in Path(manifest).read_text().splitlines() if line]
    if len(assemblies) != 464 or len({a for a, _ in assemblies}) != 464:
        raise ValueError('expected the 464-assembly HPRCv2 cohort manifest')
    rows, counts = [], {'COMPLETE': 0, 'INCOMPLETE': 0}
    for assembly, _ in sorted(assemblies):
        marker = Path(root)/'02_search'/assembly/'SUCCESS.json'
        status = 'INCOMPLETE'
        detail = 'NO_SUCCESS_MARKER'
        if marker.is_file():
            try:
                result = json.loads(marker.read_text())
                if result['status'] != 'COMPLETE' or result['monitor']['reason'] != 'completed':
                    raise ValueError('invalid completion marker')
                if result['plan_sha256'] != digest(plan):
                    raise ValueError('plan changed since completion')
                for path, spec in result['outputs'].items():
                    if digest(path) != spec['sha256']:
                        raise ValueError('output changed: '+path)
                qc = json.loads((Path(result['attempt'])/'qc-summary.json').read_text())
                if qc['status'] != 'PASS' or qc['eligible_hits'] != sum(v['eligible_hits'] for v in qc['per_guide'].values()):
                    raise ValueError('qc summary failed')
                if sum(1 for _ in open(Path(result['attempt'])/'eligible.cfd.paf')) != qc['eligible_hits']:
                    raise ValueError('PAF row count changed')
                status, detail = 'COMPLETE', 'rows=%d guides=%d' % (qc['eligible_hits'], qc['guide_count'])
            except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                status, detail = 'INCOMPLETE', str(exc)
        counts[status] += 1
        rows.append(dict(assembly=assembly, status=status, detail=detail))
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with (out.with_suffix('.tsv')).open('x') as f:
        w = csv.DictWriter(f, fieldnames=['assembly', 'status', 'detail'], delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    summary = dict(status='ALL_COMPLETE' if counts['INCOMPLETE'] == 0 else 'INCOMPLETE_ASSEMBLIES_PRESENT',
                   counts=counts, assemblies=rows,
                   caveat='Completeness of committed outputs only; not genome-wide recall.')
    out.write_text(json.dumps(summary, indent=2)+'\n')
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True)
    p.add_argument('--manifest', required=True)
    p.add_argument('--plan', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    result = verify(a.root, a.manifest, a.plan, a.out)
    print(json.dumps({k: result[k] for k in ['status', 'counts']}))
    return 0 if result['status'] == 'ALL_COMPLETE' else 1


if __name__ == '__main__':
    sys.exit(main())
