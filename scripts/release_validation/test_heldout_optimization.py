import json
from pathlib import Path
import unittest
from coverage_oracle import enumerate_sites
from resolve_bulge_discrepancies import fasta, locus

DATA=Path(__file__).resolve().parents[2]/'docs/performance_validation_a817ab5_2026-10-31_v1'
class HeldoutTests(unittest.TestCase):
    def test_independent_oracle_anchors_and_audits(self):
        for name,count in [('heldout_synthetic',223),('heldout_real',18)]:
            d=DATA/name;refs=fasta(d/'ref.fa');guides=fasta(d/'guides.fa');truth=set()
            for g,seq in guides.items():
                for contig,ref in refs.items():
                    for x in enumerate_sites(ref,seq):
                        truth.add(locus((g,contig,x['start'],x['end'],x['strand']),refs))
            actual=[]
            for line in (d/'out.paf').read_text().splitlines():
                f=line.split('\t');actual.append(locus((f[0],f[5],int(f[7]),int(f[8]),f[4]),refs))
            self.assertEqual(len(truth),count)
            self.assertEqual(truth,set(actual))
            self.assertEqual(len(actual),len(set(actual)))
            self.assertEqual(json.loads((d/'audit_summary.json').read_text())['status'],'PASS')
    def test_new_window_provenance(self):
        data=json.loads((DATA/'heldout_real/panel.json').read_text())
        self.assertEqual(len(data['windows']),4)
        self.assertTrue(all(':500001-501000' in x['region'] or ':900001-901000' in x['region'] for x in data['windows']))

if __name__=='__main__':unittest.main()
