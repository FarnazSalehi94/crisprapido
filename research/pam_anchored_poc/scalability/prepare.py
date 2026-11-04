#!/usr/bin/env python3
"""Prepare bounded, independently sourced references and bound indexes in scratch."""
import argparse,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--grch38-bgzf',required=True);p.add_argument('--samtools',required=True);p.add_argument('--builder',required=True);p.add_argument('--out',required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[3];out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
manifest=root/'research/pam_anchored_poc/compatibility/index_manifest.py'
def fasta(path):
 lines=path.read_text().splitlines();return lines[0][1:],''.join(lines[1:]).upper()
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build(name,header,seq,guide,origin):
 d=out/name;d.mkdir(exist_ok=False);ref=d/'ref.fa';ref.write_text('>'+header+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n');(d/'guides.fa').write_text('>guide\n'+guide+'\n')
 subprocess.run(['python3',str(manifest),'build','--builder',a.builder,'--reference',str(ref),'--index',str(d/'index')],check=True,stdout=(d/'build.log').open('w'),stderr=subprocess.STDOUT,timeout=180)
 return {'name':name,'source':origin,'header':header,'bases':len(seq),'Ns':seq.count('N'),'guide':guide,'reference_sha256':digest(ref),'source_sha256':origin.get('source_sha256'),'pam_gg_forward':sum(seq[i] in 'ACGT' and seq[i+1:i+3]=='GG' for i in range(len(seq)-2)),'pam_cc_reverse':sum(seq[i:i+2]=='CC' and seq[i+2] in 'ACGT' for i in range(len(seq)-2)),'index_bytes':sum(x.stat().st_size for x in d.glob('index.*'))}
items=[];guide='TCTAGCGTTTAACTTCTTCG'
for parent in ('CHM13_nonrepeat_candidate','GRCh38_repeat_candidate'):
 source=Path(a.source)/parent/'ref.fa';header,seq=fasta(source);assert len(seq)==2000000
 for size in (100000,1000000,2000000):
  # Nested prefixes preserve equal origins; headers are local to each slice.
  items.append(build(parent+'_'+str(size),header+'_prefix_'+str(size),seq[:size],guide,{'parent':header,'zero_based_offset':0,'source_sha256':digest(source),'stratum':'prior '+('repeat-enriched proxy' if 'repeat_' in parent and 'nonrepeat' not in parent else 'ordinary candidate')}))
# GRCh38 chr1 begins with a native assembly gap. Extract only this bounded region.
region='GRCh38#0#chr1:1-100000';r=subprocess.run([a.samtools,'faidx',a.grch38_bgzf,region],check=True,capture_output=True,text=True,timeout=60);lines=r.stdout.splitlines();seq=''.join(lines[1:]).upper();assert len(seq)==100000 and 'N' in seq
# Pick a resolved on-target with NGG nearest to the end of an N run.
ends=[i for i in range(1,len(seq)) if seq[i-1]=='N' and seq[i]!='N'];
choices=[(abs(i-e),i,e) for e in ends for i in range(20,len(seq)-3) if seq[i] in 'ACGT' and seq[i+1:i+3]=='GG' and set(seq[i-20:i])<=set('ACGT')]
assert choices;_,i,e=min(choices);guide=seq[i-20:i]
items.append(build('GRCh38_nativeN_100000',region,seq,guide,{'region':region,'zero_based_offset':0,'source_bgzf':a.grch38_bgzf,'n_run_end':e,'guide_pam_start':i,'guide_distance_to_N_end':i-e,'stratum':'naturally N-containing GRCh38 assembly gap edge'}))
(out/'regions.json').write_text(json.dumps(items,indent=2)+'\n');print(json.dumps(items))
