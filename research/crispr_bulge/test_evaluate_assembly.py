#!/usr/bin/env python3
"""Regression tests for saved-result auditing and deterministic sampling."""
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

from evaluate_assembly import alignment, audit, independent_requests, sample_rows


class AssemblyAuditTests(unittest.TestCase):
    def test_geometry(self):
        self.assertEqual(alignment('5=2D15=')[1:], ('DNA', 0, 2))
        self.assertEqual(alignment('5=2I13=')[1:], ('RNA', 2, 0))
        self.assertEqual(alignment('5=1I1D14=')[1:], ('mixed', 1, 1))
        with self.assertRaises(ValueError):
            alignment('20S')

    def test_preservation_and_unavailable_not_zero(self):
        raw = 'pilot\t20\t0\t20\t+\tchr\t100\t0\t22\t20\t22\t255\tnm:i:0\tcg:Z:5=2D15=\tcf:f:0.0000\n'
        with tempfile.TemporaryDirectory() as tmp:
            src, out = Path(tmp)/'in', Path(tmp)/'out'
            src.write_text(raw)
            out.write_text(raw.rstrip()+'\tbu:Z:alignment_exceeds_24_columns\n')
            rows, counts, reasons, by_geometry, _ = audit(src, out)
            self.assertEqual(counts['cfd'], 1)
            self.assertEqual(counts['crispr_bulge'], 0)
            self.assertEqual(rows[0]['cfd'], 0)
            self.assertIsNone(rows[0]['cb'])
            self.assertEqual(by_geometry['alignment_exceeds_24_columns']['DNA:I0:D2'], 1)
            out.write_text(raw.rstrip()+'\tcb:f:0\tbu:Z:alignment_exceeds_24_columns\n')
            with self.assertRaisesRegex(ValueError, 'score/status'):
                audit(src, out)
            out.write_text(raw.replace('0.0000', '0.0001').rstrip()+'\tbu:Z:scored\n')
            with self.assertRaisesRegex(ValueError, 'original fields'):
                audit(src, out)

    def test_independent_reconstruction_both_strands(self):
        sample = [dict(id=i, fields=['pilot', '20', '0', '20', strand, 'chr', '100',
                                     '0' if strand == '+' else '3', '20' if strand == '+' else '23'],
                       kind='ungapped', cigar='20=', mismatch=0)
                  for i, strand in enumerate('+-')]
        fasta = '>plus\n'+'G'*20+'AGG\n>minus\nCCT'+'C'*20+'\n'
        with tempfile.TemporaryDirectory() as tmp, patch('evaluate_assembly.subprocess.run') as run, patch(
                'evaluate_assembly.subprocess.check_output', return_value=fasta):
            work = Path(tmp)
            _, requests, strata = independent_requests(sample, work/'original.fa', work, 'G'*20)
            self.assertEqual(requests[0]['target'], requests[1]['target'])
            self.assertEqual(requests[0]['guide'], 'G'*20+'NGG')
            self.assertEqual(sum(strata.values()), 2)
            self.assertIn(str(work/'reference.fai'), run.call_args.args[0])

    def test_sampling_reproducible_includes_tails_and_all_strata(self):
        rows = [dict(id=i, cb=i/100, cfd=(100-i)/100, kind='DNA' if i % 2 else 'RNA',
                     fields=['pilot', '', '', '', '+' if i % 3 else '-'], mismatch=i % 5)
                for i in range(100)]
        sample = sample_rows(rows)
        self.assertEqual(sample, sample_rows(rows))
        ids = {r['id'] for r in sample}
        self.assertTrue(set(range(10)) <= ids)
        self.assertTrue(set(range(90, 100)) <= ids)
        key = lambda r: (r['kind'], r['fields'][4], r['mismatch'])
        self.assertEqual({key(r) for r in sample}, {key(r) for r in rows})


if __name__ == '__main__':
    unittest.main()
