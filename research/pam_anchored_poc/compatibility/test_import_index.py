#!/usr/bin/env python3
"""Integration tests against an existing small bound index; no index rebuild."""
import argparse,csv,gzip,hashlib,shutil,subprocess,sys,tempfile
from pathlib import Path
import index_manifest
from index_manifest import SUFFIXES,verify
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--index',required=True);p.add_argument('--comparator',required=True);p.add_argument('--crispr');p.add_argument('--poc');p.add_argument('--guides');p.add_argument('--expected-paf');a=p.parse_args();assert all((a.crispr,a.poc,a.guides,a.expected_paf)) or not any((a.crispr,a.poc,a.guides,a.expected_paf));here=Path(__file__).resolve().parent;root=here.parents[2];sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
with tempfile.TemporaryDirectory(prefix='crispr-import-fixture-') as t:
 d=Path(t);source=d/'source';source.mkdir();prefix=source/'idx';ref=source/'ref.fa.gz'
 with open(a.reference,'rb') as src,gzip.open(ref,'wb') as dst:shutil.copyfileobj(src,dst)
 records=[];name=None
 for line in Path(a.reference).read_text().splitlines():
  if line.startswith('>'):name=line[1:].split()[0];records.append([name,0])
  else:records[-1][1]+=len(line.strip())
 fai=Path(str(ref)+'.fai');fai.write_text(''.join('%s\t%d\t0\t80\t81\n'%(n,l) for n,l in records))
 for s in SUFFIXES:shutil.copyfile(a.index+'.'+s,str(prefix)+'.'+s)
 def write_evidence():
  with (d/'assemblies.tsv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=['assembly_id','index_prefix','fasta','fasta_sha256','fai_sha256','total_bases','contigs'],delimiter='\t');w.writeheader();w.writerow({'assembly_id':'fixture','index_prefix':str(prefix),'fasta':str(ref),'fasta_sha256':sha(ref),'fai_sha256':sha(fai),'total_bases':sum(x[1] for x in records),'contigs':len(records)})
  with (d/'components.tsv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=['assembly','path','status','sha256'],delimiter='\t');w.writeheader()
   for s in SUFFIXES:w.writerow({'assembly':'fixture','path':str(prefix)+'.'+s,'status':'READABLE','sha256':sha(str(prefix)+'.'+s)})
 def run(name):
  return subprocess.run([sys.executable,str(here/'import_index.py'),'--assembly','fixture','--assemblies',str(d/'assemblies.tsv'),'--components',str(d/'components.tsv'),'--comparator',a.comparator,'--comparator-sha256',sha(a.comparator),'--out',str(d/name)],capture_output=True,text=True,timeout=120,cwd=root)
 write_evidence();r=run('good');assert r.returncode==0,r.stderr
 assert (d/'good/index.crispr-binding.json').is_file() and not (d/'good/.incomplete').exists()
 if a.crispr:
  cmd=[a.crispr,'--search-backend','pam-anchored','--pam-search-bin',a.poc,'--columba-index',str(d/'good/index'),'-r',str(d/'good/reference.fa'),'-G',a.guides,'-t','1']
  search=subprocess.run(cmd,capture_output=True,check=True,cwd=root,timeout=120)
  assert search.stdout==Path(a.expected_paf).read_bytes(),'imported index changed integrated PAF'
 # The published binding must detect subsequent mutations of the staged files.
 for path in (d/'good/index.brt',d/'good/reference.fa'):
  path.chmod(0o600);before=path.read_bytes();path.write_bytes(before[:-1]+bytes([before[-1]^1]))
  try:verify(d/'good/reference.fa',d/'good/index')
  except ValueError:pass
  else:raise AssertionError('modified published bundle was accepted')
  path.write_bytes(before)
 verify(d/'good/reference.fa',d/'good/index')
 # Change header bytes after the first digest while preserving sequence and IDs.
 # The streamed checker must compare its own read hash against the binding.
 reference=d/'good/reference.fa';original_ref=reference.read_bytes();real_digest=index_manifest.digest;changed=[False]
 def racing_digest(path):
  value=real_digest(path)
  if Path(path)==reference and not changed[0]:
   changed[0]=True;reference.write_bytes(original_ref.replace(b'\n',b' description\n',1))
  return value
 index_manifest.digest=racing_digest
 try:
  try:verify(reference,d/'good/index')
  except ValueError as e:assert 'changed during index verification' in str(e)
  else:raise AssertionError('reference changed mid-verification was accepted')
 finally:index_manifest.digest=real_digest;reference.write_bytes(original_ref)
 verify(reference,d/'good/index')
 # A changed component cannot be made valid by matching file length.
 file=Path(str(prefix)+'.brt');original=file.read_bytes();file.write_bytes(original[:-1]+bytes([original[-1]^1]));r=run('bad_component');assert r.returncode and 'SHA256 mismatch' in r.stderr and not (d/'bad_component/index.crispr-binding.json').exists();file.write_bytes(original)
 # A new reference hash is insufficient if its content differs from the index.
 with gzip.open(ref,'rb') as f:plain=f.read()
 marker=next(i for i,x in enumerate(plain) if x in b'ACGT' and i>plain.index(b'\n'));changed=bytearray(plain);changed[marker]=ord('A') if changed[marker]!=ord('A') else ord('C')
 with gzip.open(ref,'wb') as f:f.write(changed)
 write_evidence();r=run('bad_reference');assert r.returncode and not (d/'bad_reference/index.crispr-binding.json').exists(),r.stderr
 # An incorrect trusted compressed-source hash is rejected before decompression.
 with gzip.open(ref,'wb') as f:f.write(plain)
 r=run('bad_source_digest');assert r.returncode and 'SHA256 mismatch' in r.stderr and not (d/'bad_source_digest/index.crispr-binding.json').exists(),r.stderr
 # Even a trusted digest for a truncated gzip cannot turn it into a valid input.
 raw=ref.read_bytes();ref.write_bytes(raw[:-8]);write_evidence();r=run('truncated_gzip')
 assert r.returncode and 'import failed:' in r.stderr and not (d/'truncated_gzip/index.crispr-binding.json').exists(),r.stderr
 print('import success and component/reference/source corruption rejected')
