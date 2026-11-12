#!/usr/bin/env python3
"""Bounded incident scan of INDEX_EIO assemblies for the cluster administrator.

For each assembly: full read of every index component recording the exact
first failing offset and errno; re-read of a 4 MiB window at that offset
(twice) to distinguish persistent corruption from transient errors; full
gzip+CRC verification of the source FASTA.gz (the assembly itself) so index
failures are never misreported as corrupted assemblies. Originals are only
read, never modified. One assembly per invocation.
"""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1')


def read_until_error(path, block=4*1024*1024):
    n = 0
    started = time.monotonic()
    try:
        with open(path, 'rb') as f:
            for b in iter(lambda: f.read(block), b''):
                n += len(b)
        return dict(path=path, status='READABLE', bytes_read=n, size=n, seconds=round(time.monotonic()-started, 1))
    except OSError as e:
        return dict(path=path, status='UNREADABLE', bytes_read=n, errno=e.errno, error=str(e),
                    first_failing_offset=n, seconds=round(time.monotonic()-started, 1))


def reread_window(path, offset, size=4*1024*1024):
    try:
        with open(path, 'rb') as f:
            f.seek(offset)
            data = f.read(size)
        return dict(status='READABLE', bytes=len(data))
    except OSError as e:
        return dict(status='UNREADABLE', errno=e.errno, error=str(e))


def verify_fasta(path):
    h = hashlib.sha256()
    n = 0
    started = time.monotonic()
    try:
        with gzip.open(path, 'rb') as f:
            for b in iter(lambda: f.read(4*1024*1024), b''):
                h.update(b); n += len(b)
        return dict(path=path, status='READABLE', decompressed_bytes=n, decompressed_sha256=h.hexdigest(),
                    seconds=round(time.monotonic()-started, 1))
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        return dict(path=path, status='UNREADABLE', error=str(e), decompressed_bytes=n)


def main():
    assembly = sys.argv[1]
    prefix = '/moosefs/salehi/HPRCv2_CRISPRapido/indices/{}/{}'.format(assembly, assembly)
    fasta = dict(l.split('\t') for l in
                 Path('/moosefs/salehi/HPRCv2_CRISPRapido/analysis_genomes.tsv').read_text().splitlines() if l.strip())[assembly]
    components = [read_until_error(prefix+'.'+s) for s in SUFFIXES]
    for c in components:
        if c['status'] == 'UNREADABLE':
            c['recheck_window_4mib'] = [reread_window(c['path'], c['first_failing_offset']) for _ in (1, 2)]
            c['recheck_window_zero'] = reread_window(c['path'], 0)
    source = verify_fasta(fasta)
    report = dict(assembly=assembly,
        incident_kind='persistent-unreadable-regions' if
            all(c.get('recheck_window_4mib', [{}])[0].get('status') == 'UNREADABLE' for c in components if c['status'] == 'UNREADABLE')
            else 'possibly-transient',
        failing_components=[dict(file=Path(c['path']).name, bytes_read=c['bytes_read'],
                                  first_failing_offset=c.get('first_failing_offset'), errno=c.get('errno'),
                                  error=c.get('error'), window_rechecks=c.get('recheck_window_4mib'),
                                  start_recheck=c.get('recheck_window_zero'))
                           for c in components if c['status'] == 'UNREADABLE'],
        readable_components=[dict(file=Path(c['path']).name, bytes=c['bytes_read']) for c in components if c['status'] == 'READABLE'],
        source_fasta=source,
        assembly_fastq_intact=source['status'] == 'READABLE',
        originals_modified=False)
    out = ROOT/'00_manifest/incident-evidence-v1'
    out.mkdir(exist_ok=True)
    (out/(assembly+'.json')).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(assembly=assembly, kind=report['incident_kind'],
        failing=[Path(c['path']).name for c in components if c['status'] == 'UNREADABLE'],
        fasta_ok=report['assembly_fastq_intact'])))


if __name__ == '__main__':
    main()
