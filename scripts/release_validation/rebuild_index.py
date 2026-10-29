#!/usr/bin/env python3
"""Rebuild one failed index on a NEW destination from a fully checksummed source.

Does not run CRISPRapido or queries. Input scanner evidence is mandatory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import os
import gzip


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['evidence', 'assembly', 'builder', 'release-manifest', 'output']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--prepared-status')
    p.add_argument('--time-bin', default='/home/salehi/miniconda3/bin/time')
    a = p.parse_args()
    if not Path(a.time_bin).is_file():
        raise ValueError('GNU time executable missing')
    evidence = json.loads(Path(a.evidence).read_text())
    if evidence['assembly'] != a.assembly:
        raise ValueError('wrong assembly evidence')
    failed = [r for r in evidence['components'] if r['status'] != 'READABLE' and not r['path'].endswith(('.fa.gz', '.fa.gz.fai', '.fa.gz.gzi'))]
    if not failed:
        raise ValueError('refusing unnecessary rebuild of readable index')
    source = [r for r in evidence['components'] if r['path'].endswith('.fa.gz')]
    if len(source) != 1 or source[0]['status'] != 'READABLE' or not source[0]['sha256']:
        raise ValueError('full source FASTA checksum required before rebuilding')
    out = Path(a.output)
    out.mkdir(exist_ok=False)
    for name in ['input', 'index', 'tools']:
        (out/name).mkdir()
    status = dict(status='PREPARING', assembly=a.assembly, failed_original_files=failed,
                  original_evidence=str(Path(a.evidence).resolve()), original_preserved=True)
    (out/'status.json').write_text(json.dumps(status, indent=2)+'\n')
    if a.prepared_status:
        prepared = json.loads(Path(a.prepared_status).read_text())
        old = Path(a.prepared_status).parent
        staged = old/'input'/Path(source[0]['path']).name
        plain = old/'input'/(a.assembly+'.fa')
        if sha256(staged) != source[0]['sha256'] or sha256(plain) != prepared['source_plain_sha256']:
            raise ValueError('prepared source verification failed')
        (out/'prepared_input.json').write_text(json.dumps(dict(compressed=str(staged), plain=str(plain), prior_status=a.prepared_status), indent=2)+'\n')
    else:
        staged = out/'input'/Path(source[0]['path']).name
        shutil.copyfile(source[0]['path'], str(staged))
        if sha256(staged) != source[0]['sha256']:
            raise ValueError('source changed or copy verification failed')
        plain = out/'input'/(a.assembly+'.fa')
        with gzip.open(str(staged), 'rb') as src, plain.open('xb') as dst:
            shutil.copyfileobj(src, dst, 4*1024*1024)
    builder = out/'tools/columba_build'
    shutil.copy2(a.builder, str(builder))
    builder_hash = sha256(builder)
    if builder_hash != sha256(a.builder):
        raise ValueError('builder changed during copy')
    release = json.loads(Path(a.release_manifest).read_text())
    command = [a.time_bin, '-v', '-o', str(out/'build.time'), '/usr/bin/timeout', '3600',
               release['loader'], '--library-path', release['runtime_library_path'], str(builder),
               '-f', str(plain), '-r', str(out/'index'/a.assembly)]
    status.update(status='BUILDING', command=command, source_compressed_sha256=source[0]['sha256'],
                  source_plain_sha256=sha256(plain), builder_sha256=builder_hash,
                  environment=dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    (out/'status.json').write_text(json.dumps(status, indent=2)+'\n')
    env = dict(os.environ, **status['environment'])
    with (out/'build.stdout').open('x') as stdout, (out/'build.stderr').open('x') as stderr:
        try:
            result = subprocess.run(command, env=env, stdout=stdout, stderr=stderr)
        except OSError as e:
            status.update(status='EXECUTION_FAILED', error=str(e))
            (out/'status.json').write_text(json.dumps(status, indent=2)+'\n')
            raise
    status.update(status='BUILT_NOT_YET_VERIFIED' if result.returncode == 0 else 'BUILD_FAILED',
                  exit_code=result.returncode)
    (out/'status.json').write_text(json.dumps(status, indent=2)+'\n')
    if result.returncode:
        raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
