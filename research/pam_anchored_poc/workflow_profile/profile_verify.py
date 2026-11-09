#!/usr/bin/env python3
"""Profile a full unchanged binding verification, with no guide search."""
import argparse,json,sys,time
from pathlib import Path
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'research/pam_anchored_poc/compatibility'));import index_manifest as manifest
p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--out',required=True);a=p.parse_args();stage=Path(a.stage);result={'scope':'one full existing binding verification; no search','hashes':[]};old_digest=manifest.digest;old_check=manifest.check
def digest(path):
 start=time.monotonic();v=old_digest(path);result['hashes'].append({'component':Path(path).name,'seconds':time.monotonic()-start});return v
def check(prefix,reference):
 start=time.monotonic();v=old_check(prefix,reference);result['reference_text_and_metadata_check_seconds']=time.monotonic()-start;result['reference_check']=v;return v
manifest.digest=digest;manifest.check=check;start=time.monotonic();binding=manifest.verify(stage/'reference.fa',stage/'index');result['total_seconds']=time.monotonic()-start;result['status']='passed';result['binding_reference_sha256']=binding['reference_sha256'];Path(a.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
