#!/usr/bin/env python3
"""Compare the real reconstruction expressions for supported FASTA fixtures."""
import json,subprocess,sys,tempfile
from pathlib import Path
script=Path(__file__).with_name('profile_reference_load.py')
fixtures=[b'>x\nACGTN\n>y desc\nttga\n',b'>x\r\n ACGT N \r\n\r\n>y\r\nTT\r\n',b'>x\rACGT\r>y\rNN\r',b'>x\nACGT',b'ignored\n>x\nACGT\n',b'>x\nACGT\n>x duplicate\nNN\n']
with tempfile.TemporaryDirectory() as td:
 td=Path(td)
 for i,fixture in enumerate(fixtures):
  ref=td/'ref.fa';ref.write_bytes(fixture);results=[]
  for mode in ('read_text','stream'):
   out=td/(mode+'.json');r=subprocess.run([sys.executable,str(script),'--reference',str(ref),'--mode',mode,'--out',str(out)],capture_output=True,text=True)
   results.append((r.returncode,json.loads(out.read_text()) if r.returncode==0 else None));out.unlink(missing_ok=True) if sys.version_info>=(3,8) else out.unlink() if out.exists() else None
  assert (results[0][0]==0)==(results[1][0]==0),i
  if results[0][0]==0:
   for key in ('contigs','bases','ordered_record_sha256'):assert results[0][1][key]==results[1][1][key],(i,key)
print('streaming and original reference reconstruction agree on wrapped, CRLF, CR, lowercase and duplicate-record fixtures')
