#!/usr/bin/env python3
"""Create and verify a build-time reference/index binding for vanilla indexes.

The manifest is an integrity record, not an authenticity signature. Legacy
indexes cannot prove their original N positions and must be rebuilt here.
"""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'validation'))
from index_check import check
SUFFIXES=('bwt','brt','rev.brt','sa.4','sa.bv.4','pos','sna','txt.bin','cct','meta','fsid','headerSN.bin')
def digest(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def bind_build(builder,reference,prefix):
 reference=Path(reference).resolve();prefix=Path(prefix).resolve(); manifest=Path(str(prefix)+'.crispr-binding.json')
 if manifest.exists() or any(Path(str(prefix)+'.'+s).exists() for s in SUFFIXES):raise ValueError('index prefix already exists; choose a fresh prefix')
 before=digest(reference);subprocess.run([builder,'-r',str(prefix),'-f',str(reference)],check=True)
 if digest(reference)!=before:raise ValueError('reference changed while building index')
 info=check(prefix,reference);components={s:digest(str(prefix)+'.'+s) for s in SUFFIXES}
 data={'format':'crispr-poc-index-binding-v1','reference_sha256':before,'components_sha256':components,'reference_check':info}
 with manifest.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
 return data
def verify(reference,prefix):
 prefix=Path(prefix).resolve();manifest=Path(str(prefix)+'.crispr-binding.json')
 if not manifest.is_file():raise ValueError('index binding manifest missing; import a verified existing index or build a new bound index')
 data=json.loads(manifest.read_text())
 if data.get('format')!='crispr-poc-index-binding-v1' or set(data.get('components_sha256',{}))!=set(SUFFIXES):raise ValueError('unsupported or incomplete index binding manifest')
 if digest(reference)!=data['reference_sha256']:raise ValueError('reference FASTA does not match index binding')
 for s in SUFFIXES:
  if digest(str(prefix)+'.'+s)!=data['components_sha256'][s]:raise ValueError('index component does not match binding: '+s)
 info=check(prefix,reference)
 if info['fasta_sha256']!=data['reference_sha256']:raise ValueError('reference changed during index verification')
 return data
if __name__=='__main__':
 try:
  p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
  b=sub.add_parser('build');b.add_argument('--builder',required=True);b.add_argument('--reference',required=True);b.add_argument('--index',required=True)
  v=sub.add_parser('verify');v.add_argument('--reference',required=True);v.add_argument('--index',required=True)
  a=p.parse_args();r=bind_build(a.builder,a.reference,a.index) if a.command=='build' else verify(a.reference,a.index);print(json.dumps(r,sort_keys=True))
 except (ValueError,OSError,subprocess.CalledProcessError) as e:print('index binding error: '+str(e),file=sys.stderr);sys.exit(1)
