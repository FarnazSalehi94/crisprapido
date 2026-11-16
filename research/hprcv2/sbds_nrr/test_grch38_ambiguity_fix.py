import gzip
import importlib.util
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
