import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class QueueSafetyTests(unittest.TestCase):
    def command(self, root, workers):
        return [sys.executable,str(Path(__file__).with_name('rebuild_queue.py')),
                '--root',str(root),'--evidence-dir',str(root),'--builder','/missing-builder',
                '--release-manifest','/missing-manifest','--output',str(root/'queue'),
                '--assemblies','asm','--max-builders',str(workers)]

    def test_worker_cap_is_enforced(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            result=subprocess.run(self.command(root,3),text=True,capture_output=True)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse((root/'queue').exists())

    def test_child_failure_is_not_reported_as_success(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            # A healthy set is deliberately refused by the child before any build.
            (root/'asm.json').write_text(json.dumps(dict(assembly='asm',components=[])))
            result=subprocess.run(self.command(root,1),text=True,capture_output=True)
            self.assertEqual(result.returncode,1,result.stderr)
            status=json.loads((root/'queue/summary.json').read_text())
            self.assertEqual(status['status'],'FAIL')
            self.assertNotEqual(status['finished'][0]['exit_code'],0)


if __name__=='__main__':unittest.main()
