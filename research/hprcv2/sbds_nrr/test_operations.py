import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).parent))
import campaign_snapshot as snapshot
import prepare_brunello10k as bench


class OperationsTests(unittest.TestCase):
    def test_snapshot_uses_success_selected_attempt_not_first_glob(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);names=['A'+str(i) for i in range(466)]
            (root/'analysis_genomes.tsv').write_text(''.join(n+'\tref\n' for n in names))
            r=root/'runs/sbds_ngg_nrr_v1';(r/'00_manifest').mkdir(parents=True)
            (r/'00_manifest/authoritative-manifest-v5.json').write_text(json.dumps(dict(rows=[dict(assembly=n,status='INDEX_EIO') for n in names])))
            logical=r/'02_search/A0';attempt=logical/'attempt-2';attempt.mkdir(parents=True)
            old=logical/'attempt-1';old.mkdir()
            for pam in ['NGG','NAG','NGA','NAA']:
                (old/(pam+'.summary.json')).write_text('{"audit_status":"FAIL"}')
                (attempt/(pam+'.summary.json')).write_text(json.dumps(dict(audit_status='PASS',counts=dict(rows=3,cfd_defined=3,cb_scored=2,cb_unsupported_multibase_bulge=1))))
            (logical/'SUCCESS.json').write_text(json.dumps(dict(status='COMPLETE',attempt=str(attempt))))
            s=snapshot.snapshot(root)
            self.assertEqual(s['unique_entries'],466);self.assertEqual(s['sbds_complete'],1)
            self.assertEqual(s['sbds_remaining'],465);self.assertEqual(s['totals']['hits'],12)
            self.assertEqual(s['totals']['cb_scored'],8)

    def test_benchmark_guide_selection_pin_and_uniqueness(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'guides.fa'
            def seq(i):
                chars=[]
                for _ in range(20):chars.append('ACGT'[i%4]);i//=4
                return ''.join(chars)
            p.write_text(''.join('>g'+str(i)+'\n'+seq(i)+'\n' for i in range(77441)))
            with patch.object(bench,'LIBRARY_SHA256',bench.digest(p)):
                subset=bench.select(p)
            self.assertEqual(len(subset.splitlines()),20000)
            self.assertEqual(subset.splitlines()[0],'>g0')
            self.assertEqual(subset.splitlines()[-2],'>g9999')
            with self.assertRaisesRegex(ValueError,'pin mismatch'):bench.select(p)

    def test_duplicate_cohort_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'analysis_genomes.tsv').write_text('A\tr\n'*466)
            with self.assertRaisesRegex(ValueError,'cohort identity'):snapshot.snapshot(p)


if __name__=='__main__':unittest.main()
