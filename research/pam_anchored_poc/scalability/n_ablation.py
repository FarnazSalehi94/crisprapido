#!/usr/bin/env python3
"""Create synthetic shadow using the exact resolved bases in a bound N index."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--builder',required=True);p.add_argument('--out',required=True);a=p.parse_args();source=Path(a.source);out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
root=Path(__file__).resolve().parents[3];manifest=root/'research/pam_anchored_poc/compatibility/index_manifest.py'
source_meta=next(x for x in json.loads((source/'regions.json').read_text()) if x['name']=='GRCh38_nativeN_100000');original=source/'GRCh38_nativeN_100000';orig=out/'natural';orig.mkdir()
for f in original.iterdir():
 if f.is_file():shutil.copyfile(f,orig/f.name)
lines=(original/'ref.fa').read_text().splitlines();header=lines[0];seq=''.join(lines[1:]).upper();indexed=(original/'index.txt.bin').read_bytes()[4:-1].decode('ascii');assert len(seq)==len(indexed) and all(x==y for x,y in zip(seq,indexed) if x in 'ACGT')
shadow=''.join(y if x=='N' else x for x,y in zip(seq,indexed));assert 'N' not in shadow
s=out/'shadow';s.mkdir();(s/'ref.fa').write_text(header+'\n'+'\n'.join(shadow[i:i+80] for i in range(0,len(shadow),80))+'\n');shutil.copyfile(orig/'guides.fa',s/'guides.fa')
subprocess.run(['python3',str(manifest),'build','--builder',a.builder,'--reference',str(s/'ref.fa'),'--index',str(s/'index')],stdout=(s/'build.log').open('w'),stderr=subprocess.STDOUT,check=True,timeout=180)
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest();parts=['bwt','brt','rev.brt','sa.4','sa.bv.4','pos','sna','txt.bin','cct','meta','fsid','headerSN.bin'];eq={x:sha(str(orig/'index')+'.'+x)==sha(str(s/'index')+'.'+x) for x in parts};assert all(eq.values()),eq
regions=[]
for name in ('natural','shadow'):
 m=dict(source_meta);m['name']=name;m['Ns']=seq.count('N') if name=='natural' else 0;m['reference_sha256']=sha(out/name/'ref.fa');m['source']={'kind':'native GRCh38' if name=='natural' else 'synthetic substitution of N by actual indexed bases','parent':source_meta['name']};regions.append(m)
(out/'regions.json').write_text(json.dumps(regions,indent=2)+'\n');(out/'index_component_equality.json').write_text(json.dumps(eq,indent=2)+'\n');print(json.dumps({'N_replaced':seq.count('N'),'all_components_identical':True,'components':eq}))
