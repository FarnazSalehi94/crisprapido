#!/usr/bin/env python3
"""Recover a text-only-damaged index set. Copy checksummed healthy components;
regenerate only .txt.bin from fully readable frozen source. New output only.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from verify_inputs import SUFFIXES


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for part in iter(lambda:f.read(4*1024*1024), b''):
            h.update(part)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['manifest', 'evidence', 'assembly', 'regenerate-bin', 'output']:
        p.add_argument('--'+name, required=True)
    a = p.parse_args()
    row = next(r for r in csv.DictReader(open(a.manifest),delimiter='\t') if r['assembly_id']==a.assembly)
    evidence = json.loads(Path(a.evidence).read_text())
    if evidence['assembly'] != a.assembly:
        raise ValueError('wrong evidence')
    files = {x['path']: x for x in evidence['components']}
    failed = [s for s in SUFFIXES if files[row['index_prefix']+'.'+s]['status']!='READABLE']
    if failed != ['txt.bin']:
        raise ValueError('component-only recovery requires txt.bin to be the only failed component')
    source = files[row['fasta']]
    if source['status'] != 'READABLE' or not source['sha256']:
        raise ValueError('fully readable frozen source is required')
    out = Path(a.output);out.mkdir(exist_ok=False)
    (out/'input').mkdir();(out/'index').mkdir()
    status=dict(status='PREPARING',assembly=a.assembly,method='TEXT_COMPONENT_REGENERATION',original_preserved=True,copied_components={})
    (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    staged = out/'input'/Path(row['fasta']).name
    shutil.copyfile(row['fasta'],str(staged))
    if digest(staged)!=source['sha256']:
        raise ValueError('staged source hash mismatch')
    prefix=out/'index'/Path(row['index_prefix']).name
    for suffix in SUFFIXES:
        if suffix=='txt.bin':continue
        src=row['index_prefix']+'.'+suffix;dest=str(prefix)+'.'+suffix
        shutil.copyfile(src,dest)
        value=digest(dest)
        if value!=files[src]['sha256']:
            raise ValueError('healthy component copy mismatch: '+suffix)
        status['copied_components'][suffix]=value
    command=[a.regenerate_bin,str(staged),str(prefix)+'.txt.bin',row['total_bases']]
    result=subprocess.run(command,text=True,capture_output=True)
    (out/'regeneration.stdout').write_text(result.stdout);(out/'regeneration.stderr').write_text(result.stderr)
    status.update(status='BUILT_NOT_YET_VERIFIED' if result.returncode==0 else 'REGENERATION_FAILED',command=command,source_compressed_sha256=source['sha256'],regenerate_binary_sha256=digest(a.regenerate_bin))
    (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    if result.returncode:raise SystemExit(result.returncode)


if __name__=='__main__':main()
