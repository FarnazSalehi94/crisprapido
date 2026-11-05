#!/usr/bin/env python3
"""Read-only metadata check of four panel IDs with distinct index directory names."""
import argparse,json,struct
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--index-root',required=True);p.add_argument('--panel',required=True);p.add_argument('--out',required=True);a=p.parse_args();root=Path(a.index_root);panel={x.split('\t')[0]:x.split('\t')[1] for x in Path(a.panel).read_text().splitlines()}
m={'hg002v1.1.mat_MT':'hg002v1.1.mat_MT.PanSN','hg002v1.1.pat':'hg002v1.1.pat.PanSN','HG06807_mat':'HG06807_mat_v1.0.genbank','HG06807_pat':'HG06807_pat_v1.0.genbank'};suffixes=('bwt','brt','rev.brt','sa.4','sa.bv.4','pos','sna','txt.bin','cct','meta','fsid','headerSN.bin');out=[]
for label,directory in m.items():
 prefix=root/directory/directory;item={'panel_id':label,'prefix':str(prefix)}
 try:
  item.update(meta=Path(str(prefix)+'.meta').read_text().strip().replace('\n','|'),component_sizes={s:Path(str(prefix)+'.'+s).stat().st_size for s in suffixes},binding_exists=Path(str(prefix)+'.crispr-binding.json').is_file(),status='metadata_observed')
  fa=Path(panel[label]+'.fai');expected=[(x[0],int(x[1])) for x in (line.split('\t') for line in fa.read_text().splitlines())];header=[]
  for line in Path(str(prefix)+'.headerSN.bin').read_text().splitlines():
   fields=dict(x.split(':',1) for x in line.split('\t')[1:]);header.append((fields['SN'],int(fields['LN'])))
  raw=Path(str(prefix)+'.pos').read_bytes();pos=struct.unpack('<'+'I'*(len(raw)//4),raw)
  item.update(reference=panel[label],contigs=len(expected),bases=sum(x[1] for x in expected),dictionary_matches=header==expected,offsets_match=len(pos)==len(expected)+1 and pos[0]==0 and [b-a for a,b in zip(pos,pos[1:])]==[n for _,n in expected])
 except (OSError,ValueError,KeyError,struct.error) as e:item.update(status='error',error=str(e))
 out.append(item)
Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
