#!/usr/bin/env python3
"""Positive and adversarial tests of the independent output eligibility audit."""
import subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[3];base=root/'research/pam_anchored_poc/validation/results';audit=Path(__file__).with_name('audit_paf.py')
with tempfile.TemporaryDirectory() as d:
 d=Path(d);ref=d/'ref.fa';ref.write_bytes((base/'ref.fa').read_bytes());fai=d/'ref.fa.fai'
 entries=[]
 with ref.open('rb') as f:
  while True:
   header=f.readline()
   if not header:break
   assert header.startswith(b'>');name=header[1:].strip().decode();off=f.tell();seq=f.readline();assert seq.endswith(b'\n');entries.append('%s\t%d\t%d\t%d\t%d'%(name,len(seq)-1,off,len(seq)-1,len(seq)))
 fai.write_text('\n'.join(entries)+'\n');guide='TGCCTGGAATAAATCTCACT'
 lines=[x.replace('CHM13_nonrepeat_candidate_plus_exact\t','pilot\t',1) for x in (base/'baseline.paf').read_text().splitlines() if x.startswith('CHM13_nonrepeat_candidate_plus_exact\t')];assert len(lines)==2
 paf=d/'out.paf';cmd=[sys.executable,str(audit),'--reference',str(ref),'--fai',str(fai),'--guide',guide,'--paf',str(paf)]
 def run(rows,ok):
  paf.write_text('\n'.join(rows)+'\n');p=subprocess.run(cmd,capture_output=True,text=True);assert (p.returncode==0)==ok,(p.stdout,p.stderr)
 run(lines,True);paf.write_text('\n'.join(lines)+'\n');known=cmd+['--expect-contig','CHM13_nonrepeat_candidate','--expect-strand','+','--expect-anchor','606'];assert subprocess.run(known,capture_output=True).returncode==0;assert subprocess.run(known[:-1]+['607'],capture_output=True).returncode!=0;run(lines+[lines[0]],False);run([lines[0].replace('nm:i:0','nm:i:1')],False);run([lines[0].replace('20=','19=1X')],False)
 # Modify a required PAM G in a copied reference, leaving the fixture PAF untouched.
 text=ref.read_text();text=text.replace(guide+'TGG',guide+'TNG',1);ref.write_text(text);run([lines[0]],False)
print('audit accepts eligible resolved/N sites and rejects duplicate, counts, CIGAR and invalid PAM')
