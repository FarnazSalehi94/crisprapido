#!/usr/bin/env python3
"""Independently validate selected tiny windows against original BGZF + FAI."""
import argparse,bisect,json,struct,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--cases',required=True);p.add_argument('--source-gz',required=True);p.add_argument('--samtools',required=True);p.add_argument('--stage',required=True);p.add_argument('--out',required=True);a=p.parse_args();cases=json.loads(Path(a.cases).read_text())['cases'];rows=[]
for line in Path(a.source_gz+'.fai').read_text().splitlines():
 f=line.split('\t');rows.append((f[0],int(f[1])))
posraw=(Path(a.stage)/'index.pos').read_bytes();offsets=struct.unpack('<'+'I'*(len(posraw)//4),posraw);assert [b-a for a,b in zip(offsets,offsets[1:])]==[x[1] for x in rows]
regions=[];parts=[]
for c in cases:
 i=c['contig_index'];name,L=rows[i];assert name==c['contig'] and offsets[i]+c['start']==c['global_start'];lo=c['start'];end=min(L,lo+23);r=[len(regions)];regions.append('%s:%d-%d'%(name,lo+1,end))
 if lo+23>L:
  assert i+1<len(rows);regions.append('%s:1-%d'%(rows[i+1][0],lo+23-L));r.append(len(regions)-1)
 parts.append(r)
run=subprocess.run([a.samtools,'faidx',a.source_gz]+regions,capture_output=True,text=True,timeout=120,check=True);seq=[]
for line in run.stdout.splitlines():
 if line.startswith('>'):seq.append('')
 else:seq[-1]+=line.strip().upper()
assert len(seq)==len(regions)
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];result=[]
with (Path(a.stage)/'index.txt.bin').open('rb') as text:
 for c,part in zip(cases,parts):
  original=''.join(seq[j] for j in part);assert original==c['reference_genomic'] and len(original)==23
  text.seek(4+c['global_start']);indexed=text.read(23).decode('ascii');assert indexed==c['indexed_genomic']
  assert all(x==y for x,y in zip(original,indexed) if x in 'ACGT')
  oriented=original if c['strand']=='+' else rc(original);shadow=indexed if c['strand']=='+' else rc(indexed);assert oriented==c['oriented_reference'] and shadow==c['oriented_index'] and c['guide']==shadow[:20]
  mismatch=sum(x!=y for x,y in zip(oriented[:20],shadow[:20]));assert mismatch==c['true_mismatches'];strict=oriented[20] in 'ACGT' and oriented[21:]=='GG';assert strict==c['strict_pam']
  result.append({'label':c['label'],'kind':c['kind'],'strand':c['strand'],'global_start':c['global_start'],'contig':c['contig'],'start':c['start'],'anchor':c['anchor'],'strict_pam':strict,'reference_N_target':oriented[:20].count('N'),'indexed_true_mismatches':mismatch,'source_window_agrees':True})
Path(a.out).write_text(json.dumps({'cases':result,'source_windows':len(regions),'scope':'bounded random-access source windows; not a full current source hash'},indent=2)+'\n');print('independently checked %d cases using %d source windows'%(len(result),len(regions)))
