import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from audit_emitted import (AuditError, audit, independent_cfd, load_guides,
                           matrices, parse_row, rc, validate_alignment)

ROOT = Path(__file__).resolve().parents[2]
GUIDE = 'ATCGATCGATCGATCGATCG'
SAMTOOLS = '/moosefs/guarracino/condatools/samtools/1.22.1/bin/samtools'


def record(cigar='20=', target_length=20, matches=20, mismatch=0, gaps=0, gap_size=0, cf='1.0000', strand='+', start=5):
    end = start+target_length
    return '\t'.join(map(str, ['g', 20, 0, 20, strand, 'chr', 100, start, end,
                              matches, sum(int(x) for x in __import__('re').findall(r'\d+', cigar)), 255,
                              'as:i:'+str(3*mismatch+gaps*(5+gap_size)),
                              'nm:i:'+str(mismatch), 'ng:i:'+str(gaps), 'bs:i:'+str(gap_size),
                              'cg:Z:'+cigar]+([] if cf is None else ['cf:f:'+cf])))


class AuditorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mm, cls.pm = matrices(ROOT)

    def check(self, line, target=GUIDE, pam='AGG'):
        parsed = parse_row(line, {'g': GUIDE}, {'chr': 100})
        return validate_alignment(parsed, {'g': GUIDE}, target, pam, self.mm, self.pm)

    def test_perfect_both_strands_and_all_canonical_pams(self):
        for pam in ['AGG', 'CGG', 'GGG', 'TGG']:
            for strand in ['+', '-']:
                result = self.check(record(strand=strand), GUIDE if strand=='+' else rc(GUIDE),
                                    pam if strand=='+' else rc(pam))
                self.assertEqual(result['independent_cfd'], 1)
                self.assertEqual(result['pam'], pam)

    def test_independent_mismatch_and_rna_dna_bulge_fixtures(self):
        # Stored standalone CFD fixture: position-10 mismatch score .857142857.
        target = 'ATCGATCGAACGATCGATCG'
        self.assertAlmostEqual(self.check(record('9=1X10=', matches=19, mismatch=1, cf='.8571'), target)['independent_cfd'], .857142857)
        for strand in ['+', '-']:
            for cigar, target, span in [('5=1I14=', GUIDE[:5]+GUIDE[6:], 19),
                                        ('5=1D15=', GUIDE[:5]+'A'+GUIDE[5:], 21)]:
                result = self.check(record(cigar, span, matches=19 if 'I' in cigar else 20,
                                           gaps=1, gap_size=1, cf='.4462' if 'I' in cigar else '0.0000', strand=strand),
                                    target if strand=='+' else rc(target), 'AGG' if strand=='+' else 'CCT')
                self.assertAlmostEqual(result['independent_cfd'], .446153846154 if 'I' in cigar else 0)
                self.assertEqual(result['bulges'], 1)

    def test_reject_wrong_pam_score_orientation_and_bases(self):
        for pam in ['GGA', 'ACC', 'NGG', 'GG', '']:
            with self.assertRaises(AuditError):
                self.check(record(), pam=pam)
        for line, target, pam in [(record(cf='.2593'), GUIDE, 'AGG'),
                                  (record(strand='-'), GUIDE, 'AGG'),
                                  (record(), 'T'+GUIDE[1:], 'AGG')]:
            with self.assertRaises(AuditError):
                self.check(line, target, pam)

    def test_reject_malformed_consumption_coordinates_and_tags(self):
        cases = [record().replace('cg:Z:20=', 'cg:Z:19='),
                 record().replace('cg:Z:20=', 'cg:Z:0=20='),
                 record().replace('cg:Z:20=', 'cg:Z:20M'),
                 record().replace('\t0\t20\t', '\t1\t21\t'),
                 record()+'\tcf:f:1.0', record(cf='nan'), record(cf='inf'),
                 record(cf='1.1'), record(cf=None),
                 record().replace('nm:i:0', 'nm:i:1'),
                 record().replace('as:i:0', 'as:i:9'),
                 record(strand='-', start=2), record(start=79)]
        for line in cases:
            with self.subTest(line=line):
                with self.assertRaises(AuditError):
                    self.check(line)

    def test_multiple_mismatches_and_thresholds(self):
        target = 'TTCGATCGAACGATCGATCT'
        result = self.check(record('1X8=1X9=1X', matches=17, mismatch=3, cf='.6000'), target)
        self.assertAlmostEqual(result['independent_cfd'], .5999999999)
        with self.assertRaises(AuditError):
            self.check(record('1I4=1I14=', target_length=18, matches=18, gaps=2, gap_size=1, cf='0'), GUIDE[1:5]+GUIDE[6:])

    def test_terminal_deletions_and_exact_contig_boundaries(self):
        for strand in ['+', '-']:
            for cigar, target, cf in [('1D20=', 'T'+GUIDE, '1.0000'),
                                      ('20=1D', GUIDE+'A', '.6000')]:
                result = self.check(record(cigar, 21, matches=20, gaps=1, gap_size=1,
                                           cf=cf, strand=strand),
                                    target if strand=='+' else rc(target),
                                    'AGG' if strand=='+' else 'CCT')
                self.assertAlmostEqual(result['independent_cfd'], float(cf))
        for strand, start in [('+', 0), ('-', 3)]:
            line = record(strand=strand, start=start).replace('\tchr\t100\t', '\tchr\t23\t')
            parsed = parse_row(line, {'g': GUIDE}, {'chr': 23})
            result = validate_alignment(parsed, {'g': GUIDE},
                                        GUIDE if strand=='+' else rc(GUIDE),
                                        'AGG' if strand=='+' else 'CCT', self.mm, self.pm)
            self.assertEqual(result['start0'], start)

    def test_previously_validated_real_rows(self):
        fixtures = json.loads((Path(__file__).parent/'fixtures/real_rows.json').read_text())
        self.assertTrue(any(r['paf'].split('\t')[4] == '+' for r in fixtures))
        self.assertTrue(any(r['paf'].split('\t')[4] == '-' for r in fixtures))
        for row in fixtures:
            fields = row['paf'].split('\t')
            guides = {row['guide_id']: row['guide']}
            parsed = parse_row(row['paf'], guides, {fields[5]: int(fields[6])})
            target, pam = row['oriented_target'], row['oriented_pam']
            if fields[4] == '-':
                target, pam = rc(target), rc(pam)
            result = validate_alignment(parsed, guides, target, pam, self.mm, self.pm)
            self.assertAlmostEqual(result['independent_cfd'], row['expected_cfd'], places=12)

    def test_samtools_end_to_end_batch_boundary_and_failed_output(self):
        if not Path(SAMTOOLS).exists():
            self.skipTest('explicit indexed-FASTA tool unavailable')
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            forward = 'TTTTT'+GUIDE+'AGG'+'T'*72
            reverse = 'TTTTT'+'CCT'+rc(GUIDE)+'T'*72
            (p/'ref.fa').write_text('>chr\n'+forward+'\n>rev\n'+reverse+'\n')
            subprocess.run([SAMTOOLS, 'faidx', str(p/'ref.fa')], check=True)
            (p/'guides.fa').write_text('>g\n'+GUIDE+'\n>zero\n'+'A'*20+'\n')
            negative = record(strand='-', start=8).replace('\tchr\t100\t', '\trev\t100\t')
            (p/'valid.paf').write_text(record()+'\n'+negative+'\n')
            result = audit(p/'valid.paf', p/'ref.fa', p/'guides.fa', ROOT, SAMTOOLS, p/'valid', 1)
            self.assertEqual(result['status'], 'PASS')
            self.assertEqual(result['counts']['verified'], 2)
            self.assertEqual(result['per_guide']['zero'], {})
            self.assertTrue((p/'valid/PASS').exists())
            (p/'bad.paf').write_text(record(cf='.2593')+'\n')
            result = audit(p/'bad.paf', p/'ref.fa', p/'guides.fa', ROOT, SAMTOOLS, p/'bad', 1)
            self.assertEqual(result['status'], 'FAIL')
            self.assertFalse((p/'bad/PASS').exists())
            with self.assertRaises(FileExistsError):
                audit(p/'valid.paf', p/'ref.fa', p/'guides.fa', ROOT, SAMTOOLS, p/'valid')


if __name__ == '__main__':
    unittest.main()
