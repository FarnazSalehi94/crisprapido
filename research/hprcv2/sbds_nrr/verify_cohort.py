#!/usr/bin/env python3
"""Verify completed SBDS NGG/NRR cohort assemblies and publish a report.

Reads only supervisor-completed, checksum-verified SUCCESS markers; never
interprets a missing/failed assembly as zero hits. Restarts re-run this freely.
"""
import argparse
import json
from pathlib import Path

PAMS = ('NGG', 'NAG', 'NGA', 'NAA')


def digest(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def verify(root, assemblies, out):
    report = dict(assemblies={}, missing=[], failed=[], status='PARTIAL')
    for name in assemblies:
        logical = root/'02_search'/name
        success = logical/'SUCCESS.json'
        if not success.exists():
            failed_marker = sorted(logical.glob('attempt-*/FAILED.json'))
            report['missing' if not failed_marker else 'failed'].append(
                dict(assembly=name, last_error=json.loads(failed_marker[-1].read_text())['error'] if failed_marker else 'no attempt recorded'))
            continue
        data = json.loads(success.read_text())
        for path, spec in data['outputs'].items():
            if digest(path) != spec['sha256']:
                raise ValueError('completed output changed: '+path)
        monitor_path = root/'logs'/('cohort-{}_0.json'.format(data['attempt'].rsplit('-', 1)[-1]))
        entry = dict(attempt=data['attempt'], outputs=len(data['outputs']))
        counts = {}
        for pam in PAMS:
            summary = json.loads(next(logical.glob('attempt-*/{}.summary.json'.format(pam))).read_text())
            counts[pam] = summary['counts']
        entry['counts'] = counts
        entry['totals'] = dict(hits=sum(counts[pam]['rows'] for pam in PAMS),
            cfd_defined=sum(counts[pam]['cfd_defined'] for pam in PAMS),
            cb_scored=sum(counts[pam].get('cb_scored', 0) for pam in PAMS),
            cb_na=sum(counts[pam].get(k, 0) for pam in PAMS for k in counts[pam] if k.startswith('cb_') and k != 'cb_scored'))
        report['assemblies'][name] = entry
    if not report['missing'] and not report['failed']:
        report['status'] = 'COMPLETE'
    out.write_text(json.dumps(report, indent=2)+'\n')
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--assemblies', required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(dict(status=verify(a.root, a.assemblies.split(','), a.out)['status'])))


if __name__ == '__main__':
    main()
