#!/usr/bin/env python3
"""Import a previously verified 32-bit VANILLA index without rebuilding it.

The caller must supply trusted recovery manifests and a checksum-pinned text
comparator. A new private output directory is created exclusively; failures
leave it incomplete and without a usable binding sidecar.
"""
import argparse,csv,gzip,hashlib,json,os,struct,subprocess,sys,tempfile,time
from pathlib import Path
from index_manifest import SUFFIXES,digest,verify


def copy_checked(source,dest,expected):
 h=hashlib.sha256();count=0
 with open(source,'rb') as src,open(dest,'xb') as dst:
  for block in iter(lambda:src.read(4*1024*1024),b''):
   dst.write(block);h.update(block);count+=len(block)
  dst.flush();os.fsync(dst.fileno())
 actual=h.hexdigest()
 if actual!=expected:raise ValueError('SHA256 mismatch: '+str(source))
 return {'sha256':actual,'bytes':count}

def dictionaries(prefix,fai):
 expected=[]
 for line in Path(fai).read_text().splitlines():
  f=line.split('\t');expected.append((f[0],int(f[1])))
 if not expected or len({x[0] for x in expected})!=len(expected):raise ValueError('invalid FASTA index')
 if Path(str(prefix)+'.meta').read_text()!='21\n4\nVANILLA\n':raise ValueError('unsupported index metadata')
 if Path(str(prefix)+'.fsid').read_bytes()!=struct.pack('<I',0):raise ValueError('expected one source FASTA')
 raw=Path(str(prefix)+'.pos').read_bytes()
 if len(raw)%4:raise ValueError('invalid position array')
 pos=struct.unpack('<'+'I'*(len(raw)//4),raw)
 if len(pos)!=len(expected)+1 or pos[0]!=0 or [b-a for a,b in zip(pos,pos[1:])]!=[n for _,n in expected]:raise ValueError('index offsets differ from trusted FASTA index')
 raw=Path(str(prefix)+'.sna').read_bytes();names=[];off=0
 while off<len(raw):
  if off+8>len(raw):raise ValueError('truncated sequence names')
  n=struct.unpack_from('<Q',raw,off)[0];off+=8
  if off+n>len(raw):raise ValueError('truncated sequence name')
  names.append(raw[off:off+n].decode('ascii'));off+=n
 if names!=[x[0] for x in expected]:raise ValueError('index names differ from trusted FASTA index')
 header=[]
 for line in Path(str(prefix)+'.headerSN.bin').read_text().splitlines():
  fields=dict(x.split(':',1) for x in line.split('\t')[1:]);header.append((fields['SN'],int(fields['LN'])))
 if header!=expected:raise ValueError('index header differs from trusted FASTA index')
 return {'contigs':len(expected),'bases':sum(n for _,n in expected)}

def supported_fasta(path):
 names=set();current=None;length=0
 with open(path,'rb') as f:
  for line in f:
   if line.startswith(b'>'):
    if current is not None and not length:raise ValueError('empty reference record')
    parts=line[1:].split()
    if not parts:raise ValueError('empty reference identifier')
    current=parts[0].decode('ascii')
    if current in names:raise ValueError('duplicate reference identifier')
    names.add(current);length=0
   else:
    if current is None:raise ValueError('invalid reference FASTA')
    seq=line.strip().upper()
    if any(c not in b'ACGTN' for c in seq):raise ValueError('reference contains unsupported non-ACGTN base')
    length+=len(seq)
 if current is None or not length:raise ValueError('empty reference')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--assembly',required=True);p.add_argument('--assemblies',required=True);p.add_argument('--components',required=True);p.add_argument('--comparator',required=True);p.add_argument('--comparator-sha256',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 if len(a.comparator_sha256)!=64 or digest(a.comparator)!=a.comparator_sha256:raise ValueError('comparator SHA256 mismatch')
 rows=[r for r in csv.DictReader(open(a.assemblies),delimiter='\t') if r['assembly_id']==a.assembly]
 if len(rows)!=1:raise ValueError('assembly must appear exactly once in trusted manifest')
 row=rows[0];src_prefix=row['index_prefix'];source_ref=row['fasta'];expected_ref=row['fasta_sha256'];expected_fai=row['fai_sha256']
 if len(expected_ref)!=64 or len(expected_fai)!=64:raise ValueError('trusted reference hashes missing')
 components={}
 for r in csv.DictReader(open(a.components),delimiter='\t'):
  if r['assembly']!=a.assembly:continue
  match=[s for s in SUFFIXES if r['path']==src_prefix+'.'+s]
  if len(match)!=1 or match[0] in components or r['status']!='READABLE' or len(r['sha256'])!=64:raise ValueError('invalid component evidence')
  components[match[0]]=r['sha256']
 if set(components)!=set(SUFFIXES):raise ValueError('incomplete component evidence')
 out=Path(a.out);out.mkdir(mode=0o700);start=time.monotonic();(out/'.incomplete').write_text('Import not committed; do not use.\n');prefix=out/'index';timings={};stats={}
 t=time.monotonic();stats['compressed_reference']=copy_checked(source_ref,out/'source.fa.gz',expected_ref);timings['copy_compressed_seconds']=time.monotonic()-t
 t=time.monotonic();h=hashlib.sha256()
 with gzip.open(out/'source.fa.gz','rb') as src,open(out/'reference.fa','xb') as dst:
  for block in iter(lambda:src.read(4*1024*1024),b''):dst.write(block);h.update(block)
  dst.flush();os.fsync(dst.fileno())
 ref_digest=h.hexdigest();timings['decompress_seconds']=time.monotonic()-t
 t=time.monotonic();fai=out/'reference.fa.gz.fai';stats['fai']=copy_checked(source_ref+'.fai',fai,expected_fai)
 for s in SUFFIXES:stats[s]=copy_checked(src_prefix+'.'+s,str(prefix)+'.'+s,components[s]);print('verified copy '+s,file=sys.stderr,flush=True)
 timings['copy_index_seconds']=time.monotonic()-t
 t=time.monotonic();info=dictionaries(prefix,fai);supported_fasta(out/'reference.fa')
 if info['bases']!=int(row['total_bases']) or info['contigs']!=int(row['contigs']):raise ValueError('assembly dimensions differ from trusted manifest')
 comparison=subprocess.run([a.comparator,str(out/'source.fa.gz'),str(prefix)+'.txt.bin'],capture_output=True,text=True,check=True,timeout=1800);result=json.loads(comparison.stdout)
 if result['status']!='PASS' or result['bases']!=info['bases'] or result['contigs']!=info['contigs']:raise ValueError('normalized reference/index comparison failed')
 timings['reference_text_validation_seconds']=time.monotonic()-t
 # Publish the usable binding only after all full-hash and reference checks pass.
 data={'format':'crispr-poc-index-binding-v1','reference_sha256':ref_digest,'components_sha256':components,'reference_check':{'scope':'verified import; normalized text comparison and trusted SHA256s','bases':info['bases'],'contigs':info['contigs']}}
 # Validate through the real verifier using a private sidecar namespace.
 # The consumer-visible sidecar is not published until this check succeeds.
 with tempfile.TemporaryDirectory(prefix='.verify-',dir=str(out)) as td:
  vp=Path(td)/'index'
  for s in SUFFIXES:os.symlink(str((out/('index.'+s)).resolve()),str(vp)+'.'+s)
  Path(str(vp)+'.crispr-binding.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
  t=time.monotonic();verify(out/'reference.fa',vp);timings['final_binding_verify_seconds']=time.monotonic()-t
 (out/'source.fa.gz').unlink();(out/'reference.fa.gz.fai').unlink();timings['total_seconds']=time.monotonic()-start
 summary={'assembly':a.assembly,'reference':str(out/'reference.fa'),'index':str(prefix),'source_reference':source_ref,'source_index':src_prefix,'reference_sha256':ref_digest,'source_sha256':expected_ref,'comparator_sha256':a.comparator_sha256,'components':stats,'dimensions':info,'normalized_text_comparison':result,'timings':timings,'limits':'Trusted manifests and comparator are external evidence; no FM-index search or full FM structural proof.'}
 (out/'import.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 pending=out/'index.binding.pending';pending.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n');os.replace(pending,out/'index.crispr-binding.json');(out/'.incomplete').unlink()
 for f in out.iterdir():f.chmod(0o400)
 print(json.dumps(summary,sort_keys=True))
if __name__=='__main__':
 try:main()
 except (OSError,EOFError,ValueError,KeyError,UnicodeError,struct.error,subprocess.SubprocessError) as e:print('import failed: '+str(e),file=sys.stderr);sys.exit(1)
