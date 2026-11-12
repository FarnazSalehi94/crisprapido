#!/usr/bin/env python3
"""Verify four literal PAM pilots and publish a non-overlapping NRR union."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()


def collect(root, output):
    plan_hash = digest(root/'00_manifest/pilot-plan-v1.json')
    runs = {}
    for pam in ('NGG', 'NAG', 'NGA', 'NAA'):
        s = json.loads((root/'01_pilot'/pam/'SUCCESS.json').read_text())
        if s['plan_sha256'] != plan_hash: raise ValueError('plan mismatch')
        for path, sha in s['outputs'].items():
            if digest(path) != sha: raise ValueError('output checksum mismatch: '+path)
        annotations = next(Path(p) for p in s['outputs'] if p.endswith('/annotations.tsv'))
        attempt = annotations.parents[2]
        job = attempt.name.removeprefix('attempt-')
        monitor = json.loads((root/'logs'/('monitor-'+job+'.json')).read_text())
        if monitor['reason'] != 'completed': raise ValueError('supervisor failure')
        runs[pam] = (attempt, annotations, monitor)
    output.mkdir(exist_ok=False)
    seen = set()
    counts = {}
    with (output/'nrr.both.paf').open('x') as merged, (output/'pam-classification.tsv').open('x') as labels:
        writer = csv.writer(labels, delimiter='\t')
        writer.writerow(['row', 'guide', 'contig', 'strand', 'pam_boundary', 'observed_pam', 'pam_class'])
        for pam, (attempt, annotations, monitor) in runs.items():
            n = 0
            with annotations.open() as f, (attempt/'scoring.stdout').open() as paf:
                for r in csv.DictReader(f, delimiter='\t'):
                    line = paf.readline()
                    if not line: raise ValueError('PAF shorter than annotations')
                    key = (r['guide_id'], r['contig'], r['strand'], r['end_0based'] if r['strand']=='+' else r['start_0based'])
                    if key in seen or 'N'+r['pam'][1:] != pam: raise ValueError('overlap/PAM mismatch')
                    seen.add(key); n += 1
                    merged.write(line)
                    writer.writerow([len(seen), *key, r['pam'], pam])
                if paf.readline(): raise ValueError('PAF longer than annotations')
            counts[pam] = dict(rows=n, summary=json.loads((attempt/'summary.json').read_text()), search=json.loads((attempt/'search.metrics.json').read_text()), scoring=json.loads((attempt/'scoring.metrics.json').read_text()), monitor=monitor)
    # This is a union of disjoint literal searches, not an independent combined-search recall experiment.
    report = dict(status='FOUR_LITERAL_PILOTS_AUDITED', counts=counts, total=len(seen), ngg_subset='exact copied NGG component; disjoint union checked, not independent recall evidence', neural_validation='noncanonical official computation demonstrated on four synthetic inputs only; production policy remains unsupported_pam', cohort_status='BLOCKED pending full input binding verification and noncanonical scoring review', outputs={str(p):digest(p) for p in output.iterdir() if p.is_file()})
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); collect(a.root,a.output)
