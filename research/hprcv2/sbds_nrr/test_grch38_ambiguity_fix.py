import gzip
import importlib.util
import shutil
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('fix', Path(__file__).with_name('grch38_ambiguity_fix.py'))
fix = importlib.util.module_from_spec(spec); spec.loader.exec_module(fix)


class ConversionTests(unittest.TestCase):
    def make_source(self, d):
        src = Path(d)/'source.fa.gz'
        with gzip.open(src, 'wt') as f:
            f.write('>chr1 first\nACGTRYKM\nacgtn\n>chr2 second\nNNNNBDHV\n')
        return src

    def test_derived_gzip_is_deterministic(self):
        with tempfile.TemporaryDirectory() as d:
            src = self.make_source(d)
            first, second = Path(d)/'a.gz', Path(d)/'b.gz'
            fix.convert(src, first, 'samtools')
            fix.convert(src, second, 'samtools')
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_streaming_decompression_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            src = self.make_source(d); derived = Path(d)/'d.gz'
            fix.convert(src, derived, 'samtools')
            plain = Path(d)/'d.fa'
            with gzip.open(derived, 'rb') as f, plain.open('wb') as out:
                shutil.copyfileobj(f, out, 8*1024*1024)
            with gzip.open(src, 'rt') as f:
                expected = f.read().replace('R', 'N').replace('Y', 'N').replace('K', 'N').replace('M', 'N').replace('B', 'N').replace('D', 'N').replace('H', 'N').replace('V', 'N')
            self.assertEqual(plain.read_text(), expected)

    def test_conversion_positions_and_lengths(self):
        with tempfile.TemporaryDirectory() as d:
            src = self.make_source(d); dst = Path(d)/'derived.fa.gz'
            positions, summary = fix.convert(src, dst, 'samtools')
            self.assertEqual(positions, [
                ['chr1', 4, 'R', 'N'], ['chr1', 5, 'Y', 'N'], ['chr1', 6, 'K', 'N'], ['chr1', 7, 'M', 'N'],
                ['chr2', 4, 'B', 'N'], ['chr2', 5, 'D', 'N'], ['chr2', 6, 'H', 'N'], ['chr2', 7, 'V', 'N']])
            self.assertEqual(summary['converted_positions'], 8)
            self.assertEqual(summary['headers'], 2)
            with gzip.open(dst, 'rt') as f:
                self.assertEqual(f.read(), '>chr1 first\nACGTNNNN\nacgtn\n>chr2 second\nNNNNNNNN\n')

    def test_rejects_malformed_and_unconverted(self):
        with tempfile.TemporaryDirectory() as d:
            empty = Path(d)/'empty.fa.gz'
            with gzip.open(empty, 'wt') as f:
                f.write('>only header\n')
            with self.assertRaises(ValueError):
                fix.convert(empty, Path(d)/'y.gz', 'samtools')
            clean = Path(d)/'clean.fa.gz'
            with gzip.open(clean, 'wt') as f:
                f.write('>chr1\nACGT\n')
            positions, _ = fix.convert(clean, Path(d)/'z.gz', 'samtools')
            self.assertEqual(positions, [])

    def test_rejects_invalid_character(self):
        with tempfile.TemporaryDirectory() as d:
            bad = Path(d)/'bad.fa.gz'
            with gzip.open(bad, 'wt') as f:
                f.write('>chr1\nACGTX\n')
            with self.assertRaises(ValueError):
                fix.convert(bad, Path(d)/'w.gz', 'samtools')


if __name__ == '__main__':
    unittest.main()
