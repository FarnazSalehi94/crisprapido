"""Small fixtures only; no genome/index reads. Run with python3 -m unittest discover."""
import gzip
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest


class TextComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.binary = cls.root/'compare_index_text'
        env = dict(os.environ)
        for name in ['CPATH', 'CPLUS_INCLUDE_PATH', 'C_INCLUDE_PATH', 'LIBRARY_PATH']:
            env.pop(name, None)
        subprocess.run(['/usr/bin/g++', '-O2', '-std=c++14', str(Path(__file__).with_name('compare_index_text.cpp')),
                        '-lz', '-o', str(cls.binary)], env=env, check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_fixture(self, fasta, text, declared=None):
        with tempfile.TemporaryDirectory(dir=str(self.root)) as td:
            path = Path(td)
            with gzip.open(str(path/'ref.fa.gz'), 'wb') as f:
                f.write(fasta)
            (path/'index.bin').write_bytes(struct.pack('<I', len(text) if declared is None else declared)+text)
            return subprocess.run([str(self.binary), str(path/'ref.fa.gz'), str(path/'index.bin')],
                                  text=True, capture_output=True)

    def test_multiple_records_lowercase_and_ambiguous_replacement(self):
        result = self.run_fixture(b'>first\nacgtN\n>second\nTTaa\n', b'ACGTATTAA$')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('"contigs":2', result.stdout)
        self.assertIn('"ambiguous_replacements":1', result.stdout)

    def test_chunk_crossing(self):
        seq = b'ACGT'*300000
        result = self.run_fixture(b'>first\n'+seq[:600000]+b'\n>second\n'+seq[600000:]+b'\n', seq+b'$')
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_wrong_base_fails(self):
        result = self.run_fixture(b'>first\nACGT\n', b'ACGA$')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('sequence mismatch at concatenated position 3', result.stderr)

    def test_truncated_and_wrong_length_fail(self):
        for text, declared in [(b'ACG', 5), (b'ACGT$', 99), (b'ACGTX', 5)]:
            result = self.run_fixture(b'>first\nACGT\n', text, declared)
            self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()
