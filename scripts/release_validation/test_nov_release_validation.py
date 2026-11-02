import csv
import json
from pathlib import Path
import re
import unittest
from coverage_oracle import enumerate_sites
from resolve_bulge_discrepancies import fasta,locus,overlaps

DATA=Path(__file__).resolve().parents[2]/'docs/release_validation_0aa9944_2026-11-02_v1'
class NovemberValidationTests(unittest.TestCase):
    def test_heldout_recomputed_truth_and_raw_candidate_overlap(self):
        for panel,expected in [('heldout_synthetic',212),('heldout_real',18)]:
            d=DATA/panel;refs=fasta(d/'ref.fa');guides=fasta(d/'guides.fa');truth=set();intervals=set()
            for guide,seq in guides.items():
                for contig,ref in refs.items():
                    for hit in enumerate_sites(ref,seq):
                        key=(guide,contig,hit['start'],hit['end'],hit['strand']);intervals.add(key);truth.add(locus(key,refs))
            raw=set()
            for f in ('edit.sam','exact.sam'):
                for line in (d/f).read_text().splitlines():
                    if line.startswith('@'):continue
                    x=line.split('\t');flag=int(x[1])
                    if flag&4:continue
                    s=int(x[3])-1;n=sum(int(c) for c,o in re.findall(r'(\d+)([MIDNSHP=X])',x[5]) if o in 'MDN=X')
                    raw.add((x[0],x[2],s,s+n,'-' if flag&16 else '+'))
            reported=[]
            for line in (d/'out.paf').read_text().splitlines():
                x=line.split('\t');reported.append(locus((x[0],x[5],int(x[7]),int(x[8]),x[4]),refs))
            self.assertEqual(len(truth),expected);self.assertEqual(set(reported),truth)
            self.assertEqual(len(reported),len(set(reported)))
            self.assertTrue(all(any(overlaps(x,y) for y in raw) for x in intervals))
            self.assertEqual(json.loads((d/'audit_summary.json').read_text())['status'],'PASS')
            self.assertTrue(any(a[:3]==b[:3] and a[3]<b[4] and b[3]<a[4] for a in truth for b in truth if a<b))
    def test_heldout_sequences_and_windows_are_new(self):
        synthetic=fasta(DATA/'heldout_synthetic/ref.fa')
        root=DATA.parent
        for previous in ('broader_coverage_1eb2f40_2026-10-30/synthetic/ref.fa',
                         'performance_validation_a817ab5_2026-10-31_v1/heldout_synthetic/ref.fa'):
            self.assertFalse(set(synthetic.values()) & set(fasta(root/previous).values()))
        panel=json.loads((DATA/'heldout_real/panel.json').read_text())
        self.assertTrue(all(':1300001-1301000' in x['region'] or ':1700001-1701000' in x['region'] for x in panel['windows']))

    def test_scaling_summary_and_frozen_parameters(self):
        d=DATA/'scaling';rows=list(csv.DictReader((d/'per_guide.tsv').read_text().splitlines(),delimiter='\t'))
        self.assertEqual(len(rows),16);self.assertEqual(sum(int(x['reported']) for x in rows),43)
        self.assertTrue(all(x['audit']=='PASS' and int(x['unique_raw'])>0 for x in rows))
        emitted=[json.loads(x) for x in (d/'emitted.jsonl').read_text().splitlines()]
        self.assertEqual(len(emitted),43)
        self.assertEqual(len(json.loads((d/'audit_summaries.json').read_text())),16)
        manifest=json.loads((DATA/'freeze_manifest.json').read_text())
        self.assertEqual(manifest['commit'],'0aa9944e2eae4b7f8286415ee5c3fb3667cf6097')
        self.assertEqual((manifest['search']['max_mismatches'],manifest['search']['max_bulges'],manifest['search']['max_bulge_size'],manifest['search']['min_match_fraction']),(4,1,1,.75))

if __name__=='__main__':unittest.main()
