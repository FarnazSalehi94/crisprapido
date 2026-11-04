#!/usr/bin/env python3
"""Build-time bindings reject stale references, N changes, and index tampering."""
import json,subprocess,sys
from pathlib import Path
builder,poc,out=sys.argv[1:];p=Path(out);p.mkdir(exist_ok=False);r=p/'ref.fa';g=p/'guides.fa';idx=p/'index';r.write_text('>chr\nACGTTGCACTGAGTACCGTANGG\n');g.write_text('>g\nACGTTGCACTGAGTACCGTA\n')
here=Path(__file__).resolve().parent;subprocess.run([sys.executable,str(here/'index_manifest.py'),'build','--builder',builder,'--reference',str(r),'--index',str(idx)],capture_output=True,text=True,check=True)
def run():return subprocess.run([sys.executable,str(here/'checked_search.py'),'--poc',poc,'-r',str(r),'-G',str(g),'--index',str(idx)],capture_output=True,text=True)
assert run().returncode==0
original=r.read_text();r.write_text(original.replace('N','A'));x=run();assert x.returncode!=0 and 'reference FASTA does not match' in x.stderr and not x.stdout;r.write_text(original)
f=Path(str(idx)+'.brt');data=f.read_bytes();f.write_bytes(data[:-1]+bytes([data[-1]^1]));x=run();assert x.returncode!=0 and 'index component' in x.stderr and not x.stdout;f.write_bytes(data)
assert run().returncode==0
Path(str(idx)+'.crispr-binding.json').unlink();x=run();assert x.returncode!=0 and 'manifest missing' in x.stderr and not x.stdout
# Other ambiguous reference symbols must fail, not be silently skipped.
r2=p/'other.fa';idx2=p/'otherindex';r2.write_text('>chr\nACGTTGCACTGAGTACCGTARGG\n')
subprocess.run([sys.executable,str(here/'index_manifest.py'),'build','--builder',builder,'--reference',str(r2),'--index',str(idx2)],capture_output=True,text=True,check=True)
x=subprocess.run([sys.executable,str(here/'checked_search.py'),'--poc',poc,'-r',str(r2),'-G',str(g),'--index',str(idx2)],capture_output=True,text=True)
assert x.returncode!=0 and 'reference must be nonempty A/C/G/T/N' in x.stderr and not x.stdout
print(json.dumps({'bound_build':True,'reference_N_change_rejected':True,'component_change_rejected':True,'missing_manifest_rejected':True,'unsupported_reference_symbol_rejected':True}))
