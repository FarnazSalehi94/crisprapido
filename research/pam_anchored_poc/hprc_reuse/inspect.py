#!/usr/bin/env python3
"""Read-only HPRC metadata inventory; does not read full index contents."""
import argparse,csv,json,struct,sys
from collections import Counter
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--index-root',required=True);p.add_argument('--panel',required=True);p.add_argument('--resolved',required=True);p.add_argument('--out',required=True);a=p.parse_args();root=Path(a.index_root);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
suffixes=('bwt','brt','rev.brt','sa.4','sa.bv.4','pos','sna','txt.bin','cct','meta','fsid','headerSN.bin');panel=[]
for line in Path(a.panel).read_text().splitlines():
 f=line.split('\t');assert len(f)>=2;panel.append({'assembly':f[0],'reference':f[1]})
assert len({x['assembly'] for x in panel})==len(panel)
for x in panel:
 prefix=root/x['assembly']/x['assembly'];x['prefix']=str(prefix)
 try:
  x['metadata']=Path(str(prefix)+'.meta').read_text();x['sizes']={s:Path(str(prefix)+'.'+s).stat().st_size for s in suffixes};x['binding_exists']=Path(str(prefix)+'.crispr-binding.json').is_file();x['status']='metadata_observed'
 except OSError as e:x['status']='error';x['error']=str(e)
(out/'panel_metadata.json').write_text(json.dumps(panel,indent=2)+'\n')
resolved=[]
for row in csv.DictReader(open(a.resolved),delimiter='\t'):
 prefix=row['index_prefix'];fai=Path(row['fasta']+'.fai');item={'assembly':row['assembly_id'],'index_prefix':prefix,'reference':row['fasta'],'prior_status':row['verification_status'],'prior_method':row['recovery_method']}
 try:
  data=[line.split('\t') for line in fai.read_text().splitlines()];records=[(x[0],int(x[1])) for x in data];posraw=Path(prefix+'.pos').read_bytes();pos=struct.unpack('<'+'I'*(len(posraw)//4),posraw);sna=Path(prefix+'.sna').read_bytes();off=0;names=[]
  while off<len(sna):
   n=struct.unpack_from('<Q',sna,off)[0];off+=8;names.append(sna[off:off+n].decode('ascii'));off+=n
  header=[]
  for line in Path(prefix+'.headerSN.bin').read_text().splitlines():
   f=dict(x.split(':',1) for x in line.split('\t')[1:]);header.append((f['SN'],int(f['LN'])))
  meta=Path(prefix+'.meta').read_text();txt=Path(prefix+'.txt.bin')
  with txt.open('rb') as f:
   raw=f.read(4);length=struct.unpack('<I',raw)[0];f.seek(-1,2);last=f.read(1)
  item.update(meta=meta.strip().replace('\n','|'),contigs=len(records),bases=sum(x[1] for x in records),dictionary_matches=names==[x[0] for x in records] and header==records,offsets_match=len(pos)==len(records)+1 and pos[0]==0 and [b-a for a,b in zip(pos,pos[1:])]==[n for _,n in records],text_header_matches=length==sum(x[1] for x in records)+1 and txt.stat().st_size==length+4 and last==b'$',component_sizes={s:Path(prefix+'.'+s).stat().st_size for s in suffixes},binding_exists=Path(prefix+'.crispr-binding.json').is_file(),reference_is_gzip=row['fasta'].endswith('.gz'),status='metadata_observed')
 except (OSError,ValueError,IndexError,struct.error) as e:item.update(status='error',error=str(e))
 resolved.append(item)
(out/'resolved_metadata.json').write_text(json.dumps(resolved,indent=2)+'\n')
summary={'panel_rows':len(panel),'metadata_observed':sum(x['status']=='metadata_observed' for x in panel),'metadata_values':dict(Counter(x.get('metadata','ERROR').strip().replace('\n','|') for x in panel)),'panel_binding_exists':sum(x.get('binding_exists',False) for x in panel),'panel_errors':[{'assembly':x['assembly'],'error':x.get('error')} for x in panel if x['status']!='metadata_observed'],'resolved_count':len(resolved),'resolved_observed':sum(x['status']=='metadata_observed' for x in resolved),'resolved_dictionary_matches':sum(x.get('dictionary_matches',False) for x in resolved),'resolved_offsets_match':sum(x.get('offsets_match',False) for x in resolved),'resolved_text_headers_match':sum(x.get('text_header_matches',False) for x in resolved),'resolved_binding_exists':sum(x.get('binding_exists',False) for x in resolved),'resolved_gzip':sum(x.get('reference_is_gzip',False) for x in resolved),'scope':'metadata, text length prefix and final byte only; no full component readability, hash or sequence identity check'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
