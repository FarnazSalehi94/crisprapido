#!/usr/bin/env python3
"""Cross-worker trust gate for rebuilt SBDS indexes (storage policy step 5).

Re-reads every component of a rebuilt index set from THIS worker (which must
differ from the build worker), compares against the completion manifest's
checksums, and only then marks the set TRUSTED and records the recheck.
Run on a node different from the build node; the script records both hosts.
"""
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import time

SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def main():
    assembly = sys.argv[1]
    build_host = sys.argv[2] if len(sys.argv) > 2 else None
    dest = Path('/moosefs/salehi/HPRCv2_CRISPRapido/indices_rebuilt_v1')/assembly
    manifest = json.loads((dest/'completion.json').read_text())
    recorded_host = manifest.get('build_host')
    if not recorded_host and (dest/'build-host.txt').exists():
        recorded_host = (dest/'build-host.txt').read_text().strip().split('=', 1)[1]
    if not recorded_host or (build_host and build_host != recorded_host):
        raise SystemExit('missing or inconsistent build-host evidence')
    build_host = recorded_host
    this_host = socket.gethostname()
    if this_host == build_host:
        raise SystemExit('recheck must run on a worker different from the build node: '+this_host)
    t = time.monotonic()
    for s in SUFFIXES:
        path = dest/(assembly+'.'+s)
        if sha256(path) != manifest['components_sha256'][s]:
            raise SystemExit('cross-node recheck FAILED for %s.%s' % (assembly, s))
    elapsed = time.monotonic()-t
    manifest['cross_node_recheck'] = dict(host=this_host, build_host=build_host,
        seconds=round(elapsed, 1), result='PASS', all_components_match=True)
    manifest['status'] = 'TRUSTED_CROSS_NODE'
    temporary = dest/'completion.json.recheck-tmp'
    with temporary.open('x') as f:
        json.dump(manifest, f, indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(temporary, dest/'completion.json')
    (dest/'INCOMPLETE').unlink(missing_ok=True)
    print(json.dumps(dict(assembly=assembly, status='TRUSTED_CROSS_NODE',
        recheck_host=this_host, seconds=round(elapsed, 1))))


if __name__ == '__main__':
    main()
