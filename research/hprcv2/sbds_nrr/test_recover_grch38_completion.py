import gzip
from pathlib import Path
import tempfile
import unittest
from recover_grch38_completion import verify, SUFFIXES, sha256


class RecoveryTests(unittest.TestCase):
    def test_verifies_and_rejects_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); scratch = root/'scratch'; dest = root/'dest'
            (scratch/'index').mkdir(parents=True); dest.mkdir()
            source = root/'source.gz'; source.write_bytes(gzip.compress(b'>chr1\nACGTN\n'))
            (scratch/'GRCh38.fa').write_bytes(b'>chr1\nACGTN\n')
            for s in SUFFIXES:
                for directory in (scratch/'index', dest):
                    (directory/('GRCh38.'+s)).write_bytes(s.encode())
            plain, hashes = verify(source, scratch, dest, sha256(source))
            self.assertEqual(plain, sha256(scratch/'GRCh38.fa'))
            self.assertEqual(len(hashes), len(SUFFIXES))
            with self.assertRaisesRegex(ValueError, 'source digest'):
                verify(source, scratch, dest, 'wrong')
            (dest/'GRCh38.brt').write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError, 'component mismatch'):
                verify(source, scratch, dest, sha256(source))
            (scratch/'GRCh38.fa').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'reference differs'):
                verify(source, scratch, dest, sha256(source))
