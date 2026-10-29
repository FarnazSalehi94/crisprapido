#!/usr/bin/env python3
"""Verify all original/recovered index bindings. No searches; continue on failure.

Accepts the completed sequential scan and an optional JSON override map. Newly
rebuilt components are read/checksummed fully again; readable original checksums
are reused with their explicit scan provenance. Each text comparison is full.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from verify_inputs import SUFFIXES, metadata
from scan_index_integrity import bounded_scan


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['manifest', 'integrity-summary', 'compare-bin', 'output']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--overrides')
    a = p.parse_args()
    original = {x['assembly']: x for x in json.loads(Path(a.integrity_summary).read_text())['assemblies']}
    overrides = json.loads(Path(a.overrides).read_text()) if a.overrides else {}
    rows = list(csv.DictReader(open(a.manifest), delimiter='\t'))
    names = [row['assembly_id'] for row in rows]
    if len(names) != 12 or len(set(names)) != 12 or set(names) != set(original):
        raise ValueError('expected exactly twelve distinct frozen assemblies with matching integrity evidence')
    if set(overrides) - set(names):
        raise ValueError('override contains an unknown assembly')
    out = Path(a.output)
    out.mkdir(exist_ok=False)
    summary = dict(assemblies=[], limitations='Readability, frozen-source/dictionary binding and reproducible component hashes are checked; not global search sensitivity or a proof of every FM-index data structure.')
    for row in rows:
        asm = row['assembly_id']
        entry = dict(assembly=asm, original_prefix=row['index_prefix'])
        evidence = original[asm]
        source = next(x for x in evidence['components'] if x['path'] == row['fasta'])
        entry['source_scan'] = source
        override = overrides.get(asm)
        prefix = override['index_prefix'] if override else row['index_prefix']
        entry.update(index_prefix=prefix, recovered=bool(override), components=[])
        if override:
            status = json.loads(Path(override['build_status']).read_text())
            if status['status'] != 'BUILT_NOT_YET_VERIFIED':
                entry.update(status='BLOCKED', error='builder did not finish successfully')
                summary['assemblies'].append(entry)
                continue
        if source['status'] != 'READABLE' or not source['sha256']:
            entry.update(status='BLOCKED_SOURCE', error='frozen source has no full readable checksum')
            summary['assemblies'].append(entry)
            continue
        old_by_suffix = {s: next(x for x in evidence['components'] if x['path']==row['index_prefix']+'.'+s) for s in SUFFIXES}
        for suffix in SUFFIXES:
            if override:
                result = bounded_scan(prefix+'.'+suffix, 90)
                if old_by_suffix[suffix]['status'] == 'READABLE':
                    result['matches_readable_original'] = result.get('sha256') == old_by_suffix[suffix]['sha256']
            else:
                result = dict(old_by_suffix[suffix], scan_provenance=a.integrity_summary)
            entry['components'].append(result)
        if not all(x['status']=='READABLE' for x in entry['components']):
            entry.update(status='BLOCKED_COMPONENTS', error='one or more components remain unreadable')
            summary['assemblies'].append(entry)
            continue
        new_row = dict(row, index_prefix=prefix)
        try:
            entry['metadata'] = metadata(new_row)
            fasta = override.get('source_gzip', row['fasta']) if override else row['fasta']
            if override:
                # Local staged gzip must be byte-identical to the fully scanned source.
                h = hashlib.sha256()
                with open(fasta, 'rb') as f:
                    for chunk in iter(lambda:f.read(4*1024*1024), b''):
                        h.update(chunk)
                if h.hexdigest() != source['sha256']:
                    raise ValueError('staged gzip differs from frozen source')
            command = [a.compare_bin, fasta, prefix+'.txt.bin']
            entry['command'] = command
            print(asm+' checking full source/text identity', flush=True)
            result = subprocess.run(command, text=True, capture_output=True, timeout=180)
            (out/(asm+'.text.stdout')).write_text(result.stdout)
            (out/(asm+'.text.stderr')).write_text(result.stderr)
            if result.returncode:
                raise ValueError('full indexed text comparison failed: '+result.stderr)
            entry['text_comparison'] = json.loads(result.stdout)
            if entry['text_comparison']['bases'] != int(row['total_bases']):
                raise ValueError('wrong base count')
            entry['status'] = ('REBUILT_VERIFIED' if override else 'ORIGINAL_VERIFIED')
            if any(x.get('matches_readable_original') is False for x in entry['components']):
                entry['status'] = 'REBUILT_NEEDS_REVIEW'
                entry['error'] = 'rebuilt hashes differ from readable original components'
        except (OSError, ValueError, subprocess.TimeoutExpired) as e:
            entry.update(status='BINDING_FAILED', error=str(e))
        summary['assemblies'].append(entry)
        (out/(asm+'.json')).write_text(json.dumps(entry, indent=2)+'\n')
        (out/'summary.partial.json').write_text(json.dumps(summary, indent=2)+'\n')
        print(asm+' '+entry['status'], flush=True)
    summary['status'] = 'PASS' if all(x['status'] in ['ORIGINAL_VERIFIED', 'REBUILT_VERIFIED'] for x in summary['assemblies']) else 'BLOCKED'
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    if summary['status'] == 'PASS':
        (out/'PASS').write_text('All twelve index sets passed source/text/readability verification.\n')
    return 0 if summary['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
