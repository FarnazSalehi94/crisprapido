import gzip
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest


class TextRegenerationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.binary = cls.root/'regenerate_text'
        env = dict(os.environ)
        for key in ['CPATH', 'CPLUS_INCLUDE_PATH', 'C_INCLUDE_PATH', 'LIBRARY_PATH']:
            env.pop(key, None)
        subprocess.run(['/usr/bin/g++', '-O2', '-std=c++14', str(Path(__file__).with_name('regenerate_text.cpp')), '-lz', '-o', str(cls.binary)], env=env, check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_layout_normalization_and_exclusive_creation(self):
        with tempfile.TemporaryDirectory(dir=str(self.root)) as td:
            p = Path(td)
            with gzip.open(str(p/'source.gz'), 'wb') as f:
                f.write(b'>plus\nacgtN\n>minus\nTTAA\n')
            result = subprocess.run([str(self.binary), str(p/'source.gz'), str(p/'text.bin'), '9'], capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            expected = struct.pack('<I', 10)+b'ACGTATTAA$'
            self.assertEqual((p/'text.bin').read_bytes(), expected)
            repeated = subprocess.run([str(self.binary), str(p/'source.gz'), str(p/'text.bin'), '9'], capture_output=True)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual((p/'text.bin').read_bytes(), expected)

    def test_bad_base_count_and_invalid_index_size_fail(self):
        with tempfile.TemporaryDirectory(dir=str(self.root)) as td:
            p = Path(td)
            with gzip.open(str(p/'source.gz'), 'wb') as f:
                f.write(b'>contig\nACGT\n')
            for i, n in enumerate(['5', str(2**32)]):
                result = subprocess.run([str(self.binary), str(p/'source.gz'), str(p/('bad'+str(i))), n], capture_output=True)
                self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()
