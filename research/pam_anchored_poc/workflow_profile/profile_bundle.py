#!/usr/bin/env python3
"""Time the existing snapshot-copy and unchanged verifier; never search a guide."""
import argparse,json,os,shutil,sys,tempfile,time
from pathlib import Path
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'research/pam_anchored_poc/compatibility'));import index_manifest as manifest
p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--scratch',required=True);p.add_argument('--out',required=True);a=p.parse_args();stage=Path(a.stage);stats={'scope':'copy one prebuilt bundle then execute unmodified full binding verifier three times; no search','phases':[],'verification_calls':[]};dest=Path(a.out)
def emit():dest.write_text(json.dumps(stats,indent=2)+'\n')
def timed(label,fn):
 start=time.monotonic();r=fn();stats['phases'].append({'name':label,'seconds':time.monotonic()-start});emit();return r
with tempfile.TemporaryDirectory(prefix='crispr-workflow-profile-',dir=a.scratch) as tmp:
 tmp=Path(tmp);ref=tmp/'reference.fa';idx=tmp/'index'
 def copy():
  shutil.copyfile(stage/'reference.fa',ref)
  for s in (*manifest.SUFFIXES,'crispr-binding.json'):shutil.copyfile(str(stage/'index')+'.'+s,str(idx)+'.'+s)
 timed('copy_reference_and_12_components_and_manifest',copy)
 original_digest=manifest.digest;original_check=manifest.check
 def digest(path):
  start=time.monotonic();v=original_digest(path);current['hashes'].append({'component':Path(path).name,'seconds':time.monotonic()-start});return v
 def check(prefix,reference):
  start=time.monotonic();v=original_check(prefix,reference);current['reference_text_and_metadata_check_seconds']=time.monotonic()-start;return v
 manifest.digest=digest;manifest.check=check
 for i in range(3):
  current={'hashes':[]};start=time.monotonic();manifest.verify(ref,idx);current['total_seconds']=time.monotonic()-start;stats['verification_calls'].append(current);emit()
 stats['status']='all_checks_passed';emit()
