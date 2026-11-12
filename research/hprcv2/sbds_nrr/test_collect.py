import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('nrr_collect', Path(__file__).with_name('collect.py'))
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class CollectTests(unittest.TestCase):
    def fixture(self, root):
        (root/'00_manifest').mkdir(); (root/'logs').mkdir()
        plan = root/'00_manifest/pilot-plan-v1.json'; plan.write_text('{}')
        for i, pam in enumerate(('NGG','NAG','NGA','NAA')):
            a=root/'01_pilot'/pam/('attempt-'+str(i)); (a/'neural/scoring').mkdir(parents=True)
            (a/'neural/scoring/annotations.tsv').write_text('guide_id\tcontig\tstrand\tstart_0based\tend_0based\tpam\n'+'g\tc\t+\t0\t'+str(20+i)+'\tA'+pam[1:]+'\n')
            (a/'scoring.stdout').write_text(pam+'\n')
            for f in ('summary.json','search.metrics.json','scoring.metrics.json'): (a/f).write_text('{}')
            (root/'logs'/('monitor-'+str(i)+'.json')).write_text('{"reason":"completed"}')
            (a.parent/'SUCCESS.json').write_text(json.dumps(dict(plan_sha256=c.digest(plan), outputs={str(p):c.digest(p) for p in a.rglob('*') if p.is_file()})))

    def test_union_and_distinct_classes(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t); self.fixture(r)
            report=c.collect(r,r/'union')
            self.assertEqual(report['total'],4)
            self.assertEqual((r/'union/nrr.both.paf').read_text(),'NGG\nNAG\nNGA\nNAA\n')

    def test_corruption_stops_before_publication(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t); self.fixture(r)
            (r/'01_pilot/NAG/attempt-1/scoring.stdout').write_text('changed')
            with self.assertRaisesRegex(ValueError,'checksum'): c.collect(r,r/'union')
            self.assertFalse((r/'union').exists())

    def test_monitor_failure_stops_publication(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t); self.fixture(r)
            (r/'logs/monitor-1.json').write_text('{"reason":"memory"}')
            with self.assertRaisesRegex(ValueError,'supervisor'): c.collect(r,r/'union')


if __name__ == '__main__': unittest.main()
