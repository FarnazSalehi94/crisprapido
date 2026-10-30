import json
from collections import Counter
from pathlib import Path
import unittest
from multi_exact_scanner import read_guides, scan_fasta
from resolve_bulge_discrepancies import fasta, locus

DATA = Path(__file__).resolve().parents[2]/'docs/broader_coverage_1eb2f40_2026-10-30'

class BroaderValidationTests(unittest.TestCase):
    def test_distinct_overlapping_exact_pams_on_both_strands(self):
        panel=DATA/'synthetic'
        with (panel/'ref.fa').open() as handle:
            rows=list(scan_fasta(handle,read_guides(panel/'guides.fa'),chunk_size=11))
        for strand in ('plus','minus'):
            contig='overlap_exact_'+strand
            a=[x for x in rows if x['contig']==contig and x['guide']=='overlap_a']
            b=[x for x in rows if x['contig']==contig and x['guide']=='overlap_b']
            self.assertEqual(len(a),1);self.assertEqual(len(b),1)
            self.assertEqual(a[0]['strand'], '+' if strand=='plus' else '-')
            self.assertNotEqual(a[0]['start'],b[0]['start'])
            self.assertTrue(a[0]['start'] < b[0]['end'] and b[0]['start'] < a[0]['end'])

    def test_missed_anchor_evidence_and_wfa_diagnoses(self):
        rows=json.loads((DATA/'missed_sites.json').read_text())
        self.assertEqual(len(rows),16)
        self.assertEqual(Counter(r['diagnosis'] for r in rows),{
            'all_oracle_intervals_rejected_by_selected_wfa_alignment':5,
            'accepted_intervals_start_beyond_two_base_candidate_expansion':11})
        refs={p:fasta(DATA/p/'ref.fa') for p in ('synthetic','real')}
        for row in rows:
            panel=row['panel']; target_locus=tuple(row['locus'])
            self.assertTrue(row['oracle_alignments'])
            for evidence in row['oracle_alignments']:
                key=(evidence['guide'],evidence['contig'],evidence['start'],evidence['end'],evidence['strand'])
                self.assertEqual(locus(key,refs[panel]),target_locus)
                for aln in evidence['oracle_alignments']:
                    a,b=aln['aligned_guide'],aln['aligned_target']
                    self.assertLessEqual(sum(x!=y and x!='-' and y!='-' for x,y in zip(a,b)),4)
                    self.assertLessEqual(a.count('-')+b.count('-'),1)
            if row['diagnosis'].startswith('all_'):
                self.assertTrue(all(not x['wfa_accepted'] for x in row['tested_intervals']))
            else:
                self.assertGreater(row['min_start_distance_for_wfa_accepted'],2)

    def test_reported_anchors_are_unique_and_in_oracle(self):
        for panel,expected in [('synthetic',(237,222)),('real',(20,19))]:
            d=DATA/panel;refs=fasta(d/'ref.fa');truth=json.loads((d/'truth.json').read_text())
            expected_loci={locus(tuple(x['key']),refs) for x in truth}
            actual=[]
            for line in (d/'out.paf').read_text().splitlines():
                f=line.split('\t');actual.append(locus((f[0],f[5],int(f[7]),int(f[8]),f[4]),refs))
            self.assertEqual((len(expected_loci),len(actual)),expected)
            self.assertEqual(len(actual),len(set(actual)))
            self.assertFalse(set(actual)-expected_loci)
            self.assertEqual(json.loads((d/'audit_summary.json').read_text())['status'],'PASS')

if __name__ == '__main__': unittest.main()
