import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class VerificationPanelTests(unittest.TestCase):
    def test_empty_or_incomplete_panel_cannot_claim_twelve_passed(self):
        script=Path(__file__).with_name('verify_recovery_set.py')
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            (root/'manifest.tsv').write_text('assembly_id\n')
            (root/'evidence.json').write_text(json.dumps(dict(assemblies=[])))
            result=subprocess.run([sys.executable,str(script),'--manifest',str(root/'manifest.tsv'),
                                   '--integrity-summary',str(root/'evidence.json'),
                                   '--compare-bin','/missing-helper','--output',str(root/'out')],
                                  text=True,capture_output=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('exactly twelve distinct',result.stderr)
            self.assertFalse((root/'out/PASS').exists())

    def test_blocked_sources_return_failure_exit_and_no_pass_marker(self):
        script=Path(__file__).with_name('verify_recovery_set.py')
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            rows=['assembly_id\tfasta\tindex_prefix']
            assemblies=[]
            for i in range(12):
                name='asm'+str(i);source='/not-read/'+name+'.fa.gz'
                rows.append(name+'\t'+source+'\t/not-read/'+name)
                assemblies.append(dict(assembly=name,components=[dict(path=source,status='TIMEOUT',sha256=None)]))
            (root/'manifest.tsv').write_text('\n'.join(rows)+'\n')
            (root/'evidence.json').write_text(json.dumps(dict(assemblies=assemblies)))
            result=subprocess.run([sys.executable,str(script),'--manifest',str(root/'manifest.tsv'),
                                   '--integrity-summary',str(root/'evidence.json'),
                                   '--compare-bin','/missing-helper','--output',str(root/'out')],
                                  text=True,capture_output=True)
            self.assertEqual(result.returncode,1,result.stderr)
            self.assertEqual(json.loads((root/'out/summary.json').read_text())['status'],'BLOCKED')
            self.assertFalse((root/'out/PASS').exists())


if __name__=='__main__':unittest.main()
