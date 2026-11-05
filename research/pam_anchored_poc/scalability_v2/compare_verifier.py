#!/usr/bin/env python3
"""Alternate the old and new FASTA parser with the same full binding verifier."""
import argparse,json,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--index',required=True);p.add_argument('--out',required=True);p.add_argument('--repetitions',type=int,default=5);a=p.parse_args();root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'research/pam_anchored_poc/compatibility'));import index_manifest
old=subprocess.run(['git','show','a2f2a07:research/pam_anchored_poc/validation/index_check.py'],cwd=root,text=True,capture_output=True,check=True).stdout;namespace={'__name__':'historical_index_check'};exec(compile(old,'historical_index_check.py','exec'),namespace);current=index_manifest.check;assert namespace['check'](a.index,a.reference)==current(a.index,a.reference);expected=None;results=[]
for i in range(a.repetitions):
 for kind in (('old','new') if i%2==0 else ('new','old')):
  index_manifest.check=namespace['check'] if kind=='old' else current;start=time.monotonic();info=index_manifest.verify(a.reference,a.index);elapsed=time.monotonic()-start
  if expected is None:expected=info
  assert info==expected
  results.append({'run':i,'parser':kind,'seconds':elapsed})
index_manifest.check=current;Path(a.out).write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
