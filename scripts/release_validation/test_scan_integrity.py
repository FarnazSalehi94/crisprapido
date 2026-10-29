import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from scan_index_integrity import scan_file


class IntegrityTests(unittest.TestCase):
    def test_complete_checksum_and_missing_file(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'component'
            path.write_bytes(b'abcdefgh')
            result = scan_file(path, 3)
            self.assertEqual(result['status'], 'READABLE')
            self.assertEqual(result['bytes_completely_read'], 8)
            self.assertEqual(result['sha256'], hashlib.sha256(b'abcdefgh').hexdigest())
            self.assertEqual(scan_file(Path(td)/'missing')['status'], 'UNREADABLE')

    def test_error_offset_and_no_false_full_checksum(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'component'
            path.write_bytes(b'abcdefgh')
            stream = MagicMock()
            stream.__enter__.return_value = stream
            stream.read.side_effect = [b'abcd', OSError(5, 'Input/output error')]
            with patch.object(Path, 'open', return_value=stream):
                result = scan_file(path, 4)
            self.assertEqual(result['failure_offset'], 4)
            self.assertEqual(result['requested_bytes'], 4)
            self.assertEqual(result['errno'], 5)
            self.assertEqual(result['bytes_completely_read'], 4)
            self.assertIsNone(result['sha256'])
            self.assertEqual(result['status'], 'UNREADABLE')


if __name__ == '__main__':
    unittest.main()
