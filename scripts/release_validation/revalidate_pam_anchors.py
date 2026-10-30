#!/usr/bin/env python3
"""Rebuild small fixture indexes and run frozen-parameter Columba/WFA validation."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

p = argparse.ArgumentParser(description=__doc__)
for key in ('fixtures', 'output', 'binary', 'columba', 'builder', 'loader', 'library-path'):
    p.add_argument('--'+key, required=True)
a = p.parse_args()
out = Path(a.output).resolve(); out.mkdir(exist_ok=False)
fixtures = Path(a.fixtures).resolve()
script = Path(__file__).resolve().parent
runs = []
def run(cmd, stdout, stderr):
    t = time.monotonic()
    with stdout.open('w') as o, stderr.open('w') as e:
        subprocess.run(cmd, stdout=o, stderr=e, check=True)
    return time.monotonic()-t
for name in ('indexed_synthetic', 'indexed_real_chm13', 'indexed_12_windows'):
    dest = out/name; dest.mkdir()
    for f in ('ref.fa','guides.fa'):
        shutil.copyfile(fixtures/name/f, dest/f)
    index = str(dest/'index'); ref = str(dest/'ref.fa'); guides = str(dest/'guides.fa')
    build_time = run([a.loader, '--library-path', a.library_path, a.builder, '-f', ref, '-r', index], dest/'build.log', dest/'build.stderr')
    for edit, label in ((0,'exact'), (5,'edit')):
        run([a.columba, '-r', index, '-f', guides, '-a','all','-m','edit','-e',str(edit),'-t','1','-o',str(dest/(label+'.sam'))], dest/(label+'.log'), dest/(label+'.stderr'))
    elapsed = run([a.binary,'-r',ref,'-G',guides,'-p','GG','-m','4','-b','1','-z','1','-f','0.75','-t','1','--columba-bin',a.columba,'--columba-index',index], dest/'out.paf', dest/'run.stderr')
    with (dest/'metrics.jsonl').open('w') as log:
        subprocess.run([sys.executable, str(script/'compare_indexed_coverage.py'), str(dest)], check=True, stdout=log)
    subprocess.run([sys.executable, str(script/'resolve_bulge_discrepancies.py'), str(dest), '--output', str(dest/'alignment_review.json')], check=True)
    runs.append({'panel':name, 'index_build_seconds':build_time, 'corrected_search_seconds':elapsed})
(out/'timing.json').write_text(json.dumps(runs, indent=2)+'\n')
