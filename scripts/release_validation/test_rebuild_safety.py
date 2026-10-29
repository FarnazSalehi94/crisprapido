import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).with_name('rebuild_index.py')


class RebuildSafetyTests(unittest.TestCase):
    def run_case(self, evidence, output, root):
        path = root/'evidence.json'
        path.write_text(json.dumps(evidence))
        return subprocess.run([sys.executable, str(SCRIPT), '--evidence', str(path),
                               '--assembly', 'asm', '--builder', '/missing-builder',
                               '--release-manifest', '/missing-manifest', '--output', str(output),
                               '--time-bin', sys.executable], text=True, capture_output=True)

    def test_refuse_rebuilding_healthy_index(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            result = self.run_case(dict(assembly='asm', components=[dict(path='asm.bwt', status='READABLE')]), p/'out', p)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('unnecessary rebuild', result.stderr)
            self.assertFalse((p/'out').exists())

    def test_require_fully_readable_source(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            result = self.run_case(dict(assembly='asm', components=[dict(path='asm.brt', status='UNREADABLE'), dict(path='asm.fa.gz', status='TIMEOUT', sha256=None)]), p/'out', p)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('full source FASTA checksum required', result.stderr)
            self.assertFalse((p/'out').exists())

    def test_never_overwrite_existing_output(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td); output = p/'out'; output.mkdir(); (output/'sentinel').write_text('preserve')
            result = self.run_case(dict(assembly='asm', components=[dict(path='asm.brt', status='UNREADABLE'), dict(path='asm.fa.gz', status='READABLE', sha256='required')]), output, p)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual((output/'sentinel').read_text(), 'preserve')
            self.assertEqual(len(list(output.iterdir())), 1)


if __name__ == '__main__':
    unittest.main()
