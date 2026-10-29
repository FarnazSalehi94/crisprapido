import io
import unittest
from coverage_oracle import enumerate_sites, rc
from multi_exact_scanner import scan_fasta

class MultiExactTests(unittest.TestCase):
    def test_matches_independent_window_oracle_across_line_boundaries(self):
        guides = {'a': 'ACGTACGTTGCAACGTACGA', 'b': 'AAAA', 'c': 'ACGTACGTTGCAACGTACGA'}
        seq = guides['a'] + 'AGGNN' + rc(guides['a'] + 'TGG') + 'NNAAAAAGGG'
        expected = {(name, s['start'], s['end'], s['strand'], s['pam'])
                    for name, g in guides.items() for s in enumerate_sites(seq, g, 0, False)}
        for width in (1, 3, 7, 20, 100):
            fasta = '>chr\n' + '\n'.join(seq[i:i+width] for i in range(0, len(seq), width)) + '\n'
            rows = list(scan_fasta(io.StringIO(fasta), guides))
            actual = {(x['guide'], x['start'], x['end'], x['strand'], x['pam']) for x in rows}
            self.assertEqual(actual, expected)
            self.assertEqual(len(rows), len(actual))

    def test_records_and_ambiguous_pam(self):
        rows = list(scan_fasta(io.StringIO('>a\nACGTNGG\n>b\nACGTAGG\n'), {'g': 'ACGT'}))
        self.assertEqual([(r['contig'], r['start']) for r in rows], [('b', 0)])

if __name__ == '__main__': unittest.main()
