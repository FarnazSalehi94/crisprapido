"""Tests for the bounded workflow validator, not another search implementation."""
from types import SimpleNamespace
import unittest

from validate_workflow import check_oracle, independent_cfd_check, load_tables, parse_rows


class WorkflowAuditTests(unittest.TestCase):
    raw = 'g\t20\t0\t20\t+\tchr\t23\t0\t20\t20\t20\t255\tnm:i:0\tcg:Z:20=\tcf:f:1.0000\n'

    def test_duplicates_rejected_but_query_ids_distinct(self):
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            parse_rows(self.raw*2)
        self.assertEqual(len(parse_rows(self.raw+self.raw.replace('g\t', 'other\t', 1))), 2)

    def test_missing_unexpected_and_changed_representative_rejected(self):
        expected = {('g', 'chr', '+', 20): SimpleNamespace(start=0, end=20, cigar='20=', mm=0, score=1)}
        self.assertEqual(check_oracle(self.raw, expected), 1)
        for wrong in ('', self.raw.replace('chr', 'wrong'), self.raw.replace('20=', '19=1X'),
                      self.raw.replace('1.0000', '0.0000')):
            with self.assertRaises(ValueError):
                check_oracle(wrong, expected)

    def test_unavailable_not_promoted_to_zero(self):
        absent = self.raw.replace('\tcf:f:1.0000', '')
        expected = {('g', 'chr', '+', 20): SimpleNamespace(start=0, end=20, cigar='20=', mm=0, score=None)}
        self.assertEqual(check_oracle(absent, expected), 1)
        with self.assertRaises(ValueError):
            check_oracle(self.raw.replace('1.0000', '0.0000'), expected)

    def test_independent_cfd_replay_detects_score_mutation(self):
        guide = 'GAAACAGTCGATTTTATCAC'
        self.assertEqual(independent_cfd_check(self.raw, {'g': guide}, {'chr': guide+'AGG'}, load_tables()), {'defined': 1})
        with self.assertRaisesRegex(ValueError, 'CFD'):
            independent_cfd_check(self.raw.replace('1.0000', '0.9000'), {'g': guide}, {'chr': guide+'AGG'}, load_tables())


if __name__ == '__main__':
    unittest.main()
