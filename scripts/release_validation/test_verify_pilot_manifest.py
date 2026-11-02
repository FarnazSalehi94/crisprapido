import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from verify_pilot_manifest import SUFFIXES

SCRIPT=Path(__file__).with_name('verify_pilot_manifest.py')
class PilotManifestTests(unittest.TestCase):
    def test_all_components_and_tamper_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);prefix=str(d/'index');fasta=d/'ref.fa';binary=d/'bin'
            def put(path,content):
                Path(path).write_bytes(content);return hashlib.sha256(content).hexdigest()
            bsha=put(binary,b'binary');fsha=put(fasta,b'>chr\nACGT\n');fai=put(str(fasta)+'.fai',b'index');gzi=put(str(fasta)+'.gzi',b'offset')
            assembly={'assembly_id':'test','index_prefix':prefix,'fasta':str(fasta),'fasta_bytes':str(fasta.stat().st_size),'fasta_sha256':fsha,'fai_sha256':fai,'gzi_sha256':gzi,'verification_status':'REBUILT_VERIFIED'}
            rows=[]
            for suffix in SUFFIXES:
                path=prefix+'.'+suffix;h=put(path,suffix.encode());rows.append({'assembly':'test','path':path,'sha256':h,'size':str(Path(path).stat().st_size),'status':'READABLE'})
            for name,values in [('assemblies.tsv',[assembly]),('checks.tsv',rows)]:
                with (d/name).open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(values[0]),delimiter='\t');w.writeheader();w.writerows(values)
            (d/'freeze.json').write_text(json.dumps({'binary_sha256':bsha}))
            base=[sys.executable,str(SCRIPT),'--assemblies',str(d/'assemblies.tsv'),'--checksums',str(d/'checks.tsv'),'--binary',str(binary),'--freeze',str(d/'freeze.json'),'--rate-mib','0','--expected-assemblies','1']
            first=subprocess.run(base+['--output',str(d/'pass')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            self.assertEqual(first.returncode,0,first.stderr);self.assertTrue((d/'pass/PASS').exists())
            Path(prefix+'.brt').write_bytes(b'bad')
            failed=subprocess.run(base+['--output',str(d/'fail')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            self.assertNotEqual(failed.returncode,0);self.assertFalse((d/'fail/PASS').exists())
            self.assertEqual(json.loads((d/'fail/summary.json').read_text())['status'],'FAIL')
if __name__=='__main__':unittest.main()
