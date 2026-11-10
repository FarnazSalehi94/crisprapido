#!/usr/bin/env python3
"""Explicit setup only: sparse, pinned upstream checkout into a NEW directory.

No network requests occur during scoring. No datasets.zip, submodules, or genomic
reference files are downloaded by this helper.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess

from adapter import digest
from worker import verify_models

p=argparse.ArgumentParser()
p.add_argument('--out',type=Path,required=True)
p.add_argument('--ca-file',type=Path,help='optional trusted CA bundle; never disable certificate checks')
a=p.parse_args()
if a.out.exists():raise SystemExit('output already exists; choose a new directory')
spec=json.loads(Path(__file__).with_name('model_manifest.json').read_text())
env={**os.environ,'GIT_LFS_SKIP_SMUDGE':'1'}
base=['git']
if a.ca_file:base+=['-c','http.sslCAInfo='+str(a.ca_file.resolve())]
def run(args):subprocess.run(base+list(map(str,args)),env=env,check=True)
run(['clone','--depth','1','--filter=blob:none','--no-checkout',spec['repository'],a.out])
if a.ca_file:run(['-C',a.out,'config','http.sslCAInfo',a.ca_file.resolve()])
run(['-C',a.out,'fetch','--depth','1','origin',spec['commit']])
paths=['OT_deep_score_src','train_and_predict_scripts','models']
paths+=sorted({str(Path(m['path']).parent) for m in spec['models']})
# Start with root files only, then detach at the pinned tree before choosing
# paths; do not depend on the remote default branch retaining the same layout.
run(['-C',a.out,'sparse-checkout','init','--cone'])
run(['-C',a.out,'checkout','--detach',spec['commit']])
run(['-C',a.out,'sparse-checkout','set','--cone',*paths])
for path,value in spec['source_sha256'].items():
    if digest(a.out/path)!=value:raise ValueError('upstream source digest mismatch: '+path)
if digest(a.out/'LICENSE')!=spec['license_sha256']:raise ValueError('upstream license digest mismatch')
verify_models(a.out,spec)
for m in spec['models']:
    if digest(a.out/m['pickle_path'])!=m['pickle_sha256']:raise ValueError('upstream pickle digest mismatch')
print('Verified pinned source and five model components; weights remain outside CRISPRapido.')
