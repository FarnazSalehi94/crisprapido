#!/usr/bin/env python3
"""Extract two 3 Mb windows using indexed FASTA; build bound indexes under limits."""
import argparse,hashlib,json,os,signal,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--chm13',required=True);p.add_argument('--grch38',required=True);p.add_argument('--samtools',required=True);p.add_argument('--builder',required=True);p.add_argument('--out',required=True);a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[3];manifest=root/'research/pam_anchored_poc/compatibility/index_manifest.py';records=[]
for assembly,genome,reg in [('CHM13',a.chm13,'CHM13#0#chr1:100000001-103000000'),('GRCh38',a.grch38,'GRCh38#0#chr9:34000001-37000000')]:
 d=out/assembly;d.mkdir();r=subprocess.run([a.samtools,'faidx',genome,reg],capture_output=True,text=True,check=True,timeout=60);lines=r.stdout.splitlines();seq=''.join(lines[1:]).upper();assert len(seq)==3000000 and set(seq)<=set('ACGTN')
 ref=d/'ref.fa';ref.write_text('>'+reg+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n');(d/'guides.fa').write_text('>guide\nTCTAGCGTTTAACTTCTTCG\n')
 start=time.monotonic();cmd=['/bin/bash','-c','ulimit -v 2097152; exec "$@"','bounded','python3',str(manifest),'build','--builder',a.builder,'--reference',str(ref),'--index',str(d/'index')]
 with (d/'build.log').open('w') as log:
  proc=subprocess.Popen(cmd,stdout=log,stderr=log,start_new_session=True)
  try:code=proc.wait(timeout=180)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait();raise RuntimeError('index build exceeded 180s')
  if code:raise RuntimeError('index build failed: '+assembly)
 build=time.monotonic()-start;start=time.monotonic();subprocess.run(['python3',str(manifest),'verify','--reference',str(ref),'--index',str(d/'index')],stdout=subprocess.DEVNULL,check=True,timeout=60);verify=time.monotonic()-start
 records.append({'name':assembly,'region':reg,'bases':len(seq),'Ns':seq.count('N'),'ref_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'index_bytes':sum(x.stat().st_size for x in d.glob('index.*')),'index_build_and_binding_seconds':build,'independent_verify_seconds':verify,'guide':'TCTAGCGTTTAACTTCTTCG'})
(out/'regions.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records))
