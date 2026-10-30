import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from multi_exact_scanner import read_guides, scan_fasta

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'docs/coverage_oracles_2026-10-29_v1'
SCRIPT = Path(__file__).with_name('compare_indexed_coverage.py')

class IndexedComparisonTests(unittest.TestCase):
    def test_recorded_fixture_recomputes_and_exact_scanner_agrees(self):
        for fixture, counts in [('indexed_synthetic', (4, 57)), ('indexed_real_chm13', (2, 7)), ('indexed_12_windows', (1, 8))]:
            with self.subTest(fixture=fixture), tempfile.TemporaryDirectory() as tmp:
                source = DATA / fixture
                for name in ('ref.fa', 'guides.fa', 'exact.sam', 'edit.sam', 'out.paf'):
                    shutil.copyfile(source / name, Path(tmp) / name)
                result = subprocess.run([sys.executable, str(SCRIPT), tmp], check=True,
                                        stdout=subprocess.PIPE, universal_newlines=True)
                metrics = [json.loads(line) for line in result.stdout.splitlines()]
                self.assertEqual(metrics[-2]['truth'], counts[0])
                self.assertEqual(metrics[-2]['found'], counts[0])
                self.assertEqual(metrics[-2]['unexpected_strict'], [])
                self.assertEqual(metrics[-1]['truth_loci'], counts[1])
                with (source / 'ref.fa').open() as handle:
                    exact = list(scan_fasta(handle, read_guides(source / 'guides.fa')))
                truth = json.loads((Path(tmp) / 'truth.json').read_text())
                expected = {tuple(x['key']) for x in truth if 'exact' in x['categories']}
                observed = {(x['guide'], x['contig'], x['start'], x['end'], x['strand']) for x in exact}
                self.assertEqual(observed, expected)
                self.assertEqual(len(exact), len(observed))
                discrepancies = json.loads((Path(tmp) / 'discrepancies.json').read_text())
                self.assertTrue(all(x['overlap_raw_edit'] and x['overlap_final'] for x in discrepancies))

if __name__ == '__main__':
    unittest.main()
