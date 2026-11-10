"""Pinned offset reuse must never bypass reference/index identity checks."""
import json
from pathlib import Path
import tempfile
import unittest

from adapter import IndexedFasta, digest
from reference_index import PinnedIndexedFasta


class PinnedIndexTests(unittest.TestCase):
    def test_valid_reuse_and_wrong_or_changed_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); ref = root/'ref.fa'; index = root/'index.json'
            ref.write_bytes(b'>one\r\nACGT\r\nTN\r\n>two\nACGT')
            original = IndexedFasta(ref)
            sha = original.sha256
            index.write_text(json.dumps(dict(format='crisprapido-fasta-offsets-v1',
                reference_sha256=sha, reference_bytes=ref.stat().st_size, rows=original.rows)))
            pin = digest(index)
            cached = PinnedIndexedFasta(ref, index, pin, sha)
            for name, (length, *_) in original.rows.items():
                for start in range(length):
                    self.assertEqual(cached.fetch(name, start, length), original.fetch(name, start, length))
            cached.close(); original.close()
            for pins in ((None, sha), (pin, None), ('0'*64, sha), (pin, '0'*64)):
                with self.assertRaises(ValueError):
                    PinnedIndexedFasta(ref, index, *pins)
            ref.write_bytes(ref.read_bytes().replace(b'ACGT', b'TGCA'))
            with self.assertRaisesRegex(ValueError, 'reference checksum'):
                PinnedIndexedFasta(ref, index, pin, sha)
            index.write_text(index.read_text()+' ')
            with self.assertRaisesRegex(ValueError, 'index checksum'):
                PinnedIndexedFasta(ref, index, pin, sha)


if __name__ == '__main__':
    unittest.main()
