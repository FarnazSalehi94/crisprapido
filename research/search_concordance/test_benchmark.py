"""Audit/normalization regressions; no search engine is mocked into concordance."""
from types import SimpleNamespace
import unittest

from benchmark import load_tables, rc, rows, site_map, witness_at, witness_sam

G = 'GAAACAGTCGATTTTATCAC'


def line(cigar='20=', hi=20, length=23, strand='+', lo=0, qlo=0, qhi=20, legacy=False):
    fields = ['Guide', '20', str(qlo), str(qhi), strand, 'chr', str(length), str(lo), str(hi), '20']
    if not legacy:
        fields.append('20')
    fields += ['255', 'nm:i:0', 'ng:i:0', 'cg:Z:'+cigar, 'cf:f:1.0000']
    return '\t'.join(fields)+'\n'


class ConcordanceAuditTests(unittest.TestCase):
    def setUp(self):
        self.tables = load_tables()

    def test_anchor_is_pam_adjacent_boundary_not_overlapping_interval(self):
        plus = rows(line(), G, {'chr': G+'AGG'}, 3, 1, self.tables)[0]
        minus = rows(line(hi=23, lo=3, strand='-'), G, {'chr': rc(G+'AGG')}, 3, 1, self.tables)[0]
        self.assertEqual(plus['key'], ['chr', '+', 20])
        self.assertEqual(minus['key'], ['chr', '-', 3])
        self.assertTrue(plus['comparable_record'] and minus['comparable_record'])
        self.assertTrue(plus['cfd_agrees'] and minus['cfd_agrees'])

    def test_legacy_schema_is_not_silently_repaired_to_paf(self):
        record = rows(line(legacy=True), G, {'chr': G+'AGG'}, 3, 1, self.tables)[0]
        self.assertTrue(record['eligible'])
        self.assertFalse(record['comparable_record'])
        self.assertEqual(record['schema'], 'legacy11_missing_block_length')
        self.assertIn('missing_mandatory_block_length', record['export_errors'])

    def test_wrong_gap_span_and_query_coordinates_are_reported_not_rewritten(self):
        seq = G+'T'+'AGG'
        bad = rows(line(cigar='20=1D', hi=20, length=24), G, {'chr': seq}, 3, 1, self.tables)[0]
        self.assertEqual(bad['key'], ['chr', '+', 20])
        self.assertIn('cigar_reference_span_disagrees_with_export', bad['alignment_errors'])
        original = rows(line(cigar='20=1D', hi=21, length=24, qlo=1, qhi=21), G, {'chr': seq}, 3, 1, self.tables)[0]
        self.assertTrue(original['eligible'])
        self.assertFalse(original['comparable_record'])
        self.assertIn('query_coordinates_not_full_guide', original['export_errors'])

    def test_mixed_witness_and_unavailable_score(self):
        seq = G[:5]+'C'+G[6:]+'AGG'
        key = ('chr', '+', 20)
        self.assertFalse(witness_at({'chr': seq}, G, key, 0, 1, self.tables))
        witnesses = witness_at({'chr': seq}, G, key, 0, 2, self.tables)
        self.assertTrue(witnesses)
        self.assertTrue(all(w.score is None for w in witnesses))
        self.assertFalse(witness_at({'chr': seq}, G, ('chr', '+', 1000), 0, 2, self.tables))

    def test_oracle_sam_has_required_tags_and_reference_oriented_negative_cigar(self):
        hit = SimpleNamespace(start=7, cigar='2=1D18=', mm=0, ins=0, dele=1)
        positive = witness_sam(('chr', '+', 28), hit)
        negative = witness_sam(('chr', '-', 7), hit)
        self.assertIn('\t2=1D18=\t', positive)
        self.assertIn('\t18=1D2=\t', negative)
        self.assertTrue(positive.endswith('\tNM:i:1\tAS:i:1\n'))
        self.assertEqual(negative.split('\t')[1], '16')

    def test_duplicate_records_are_preserved_for_audit(self):
        records = rows(line()*2, G, {'chr': G+'AGG'}, 3, 1, self.tables)
        self.assertEqual(len(site_map(records)[('chr', '+', 20)]), 2)


if __name__ == '__main__':
    unittest.main()
