import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from resolve_bulge_discrepancies import analyze, locus

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'docs/coverage_oracles_2026-10-29_v1'
TOOLS = Path(__file__).resolve().parent

class LocusReviewTests(unittest.TestCase):
    def test_distinct_pams_are_distinct_loci_even_with_overlap(self):
        refs = {'chr': 'A'*30 + 'AGGG' + 'A'*30}
        a = ('g','chr',10,30,'+')
        b = ('g','chr',11,31,'+')
        self.assertNotEqual(locus(a, refs), locus(b, refs))
        self.assertEqual(locus(a, refs)[5], 'AGG')
        self.assertEqual(locus(b, refs)[5], 'GGG')

    def test_review_recomputes_without_dropping_discrepancies(self):
        cases = [('indexed_synthetic', 46, 40, 6, 33, 28),
                 ('indexed_real_chm13', 10, 8, 2, 3, 2),
                 ('indexed_12_windows', 9, 6, 3, 4, 2)]
        for name, total, same, distinct, expected, recovered in cases:
            with self.subTest(name=name):
                rows, summary = analyze(DATA/name)
                self.assertEqual(len(rows), total)
                self.assertEqual(summary['class_counts'], {'B':same, 'C':distinct})
                self.assertEqual((summary['oracle_loci'], summary['recovered_loci']), (expected, recovered))
                self.assertEqual(summary['raw_replay_errors'], 0)
                self.assertTrue(all(x['oracle_alignments'] and x['overlapping_raw'] and x['overlapping_final'] for x in rows))
                self.assertTrue(all(x['same_pam_reported'] for x in rows if x['classification']=='B'))
                self.assertTrue(all(not x['same_pam_reported'] for x in rows if x['classification']=='C'))

    def test_genome_exact_comparison_on_synthetic_indexed_fixture(self):
        samtools = Path('/export/local/home/salehi/verkko_env/bin/samtools')
        if not samtools.exists(): self.skipTest('samtools unavailable')
        fixture = DATA/'indexed_synthetic'
        with tempfile.TemporaryDirectory() as tmp:
            ref = Path(tmp)/'ref.fa'; ref.write_text((fixture/'ref.fa').read_text())
            subprocess.run([str(samtools), 'faidx', str(ref)], check=True)
            oracle = Path(tmp)/'oracle.jsonl'; output = Path(tmp)/'out.json'
            with oracle.open('w') as out:
                subprocess.run([sys.executable, str(TOOLS/'multi_exact_scanner.py'), '--fasta', str(ref),
                    '--guides', str(fixture/'guides.fa')], stdout=out, check=True)
            subprocess.run([sys.executable, str(TOOLS/'compare_genome_exact.py'), '--fasta', str(ref),
                '--guides', str(fixture/'guides.fa'), '--oracle', str(oracle), '--sam', str(fixture/'exact.sam'),
                '--samtools', str(samtools), '--output', str(output)], stdout=subprocess.PIPE, check=True)
            result = json.loads(output.read_text())
            self.assertEqual((result['oracle_loci'],result['recovered']), (4,4))
            self.assertEqual(result['missed'], [])
            self.assertEqual(result['unexpected'], [])

if __name__ == '__main__': unittest.main()
