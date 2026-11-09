#!/usr/bin/env python3
"""Check independent multiguide auditing on a small existing bounded panel."""
import json,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[3];data=root/'research/pam_anchored_poc/validation/results';script=Path(__file__).with_name('audit_batch.py');lines=(data/'guides.fa').read_text().splitlines();guides={lines[i][1:]:lines[i+1] for i in range(0,len(lines),2)};found={};rows=[]
for line in (data/'baseline.paf').read_text().splitlines():
 f=line.split('\t');q=f[0]
 if q not in found and 'cg:Z:20=' in f:
  found[q]={'contig':f[5],'strand':f[4],'anchor':int(f[8]) if f[4]=='+' else int(f[7])-3};rows.append(line)
 if len(found)==3:break
assert len(found)==3
with tempfile.TemporaryDirectory() as tmp:
 tmp=Path(tmp);ref=tmp/'ref.fa';ref.write_bytes((data/'ref.fa').read_bytes());fai=tmp/'ref.fai';entries=[]
 with ref.open('rb') as f:
  while True:
   header=f.readline()
   if not header:break
   assert header.startswith(b'>');name=header[1:].strip().decode();offset=f.tell();seq=f.readline();assert seq.endswith(b'\n');entries.append('%s\t%d\t%d\t%d\t%d'%(name,len(seq)-1,offset,len(seq)-1,len(seq)))
 fai.write_text('\n'.join(entries)+'\n');gf=tmp/'guides.fa';gf.write_text(''.join('>'+q+'\n'+guides[q]+'\n' for q in found));targets=tmp/'targets.json';targets.write_text(json.dumps(found));paf=tmp/'out.paf';cmd=[sys.executable,str(script),'--reference',str(ref),'--fai',str(fai),'--guides',str(gf),'--targets',str(targets),'--paf',str(paf)]
 def run(lines,success):
  paf.write_text('\n'.join(lines)+'\n');r=subprocess.run(cmd,capture_output=True,text=True);assert (r.returncode==0)==success,(r.stdout,r.stderr)
 run(rows,True);run(rows+[rows[0]],False);run(rows[:-1],False);run([rows[0].replace('nm:i:0','nm:i:1')]+rows[1:],False)
print('three-guide audit passed; duplicate, missing exact target and altered mismatch count rejected')
