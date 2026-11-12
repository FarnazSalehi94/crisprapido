import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, str(HERE/(name+'.py')))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


worker = load('bench16_worker')
analyze = load('analyze16')


class PlanTests(unittest.TestCase):
    def plan(self, tmp='/tmp'):
        assemblies = {}
        for name, ref in [('HG00735_pat', 'a'*64), ('HG01496_pat', 'b'*64)]:
            root = Path(tmp)/name
            root.mkdir(parents=True, exist_ok=True)
            (root/'binding.json').write_text(json.dumps(dict(reference_sha256=ref)))
            assemblies[name] = dict(reference=str(root/'ref.fa'), reference_sha256=ref,
                                    index=str(root/'index'), binding=str(root/'binding.json'),
                                    fai=str(root/'ref.fa.fai'), scoring_reference_index=str(root/'scoring-index.json'))
        pins = {v: '0'*64 for spec in assemblies.values() for v in spec.values()}
        for spec in assemblies.values():
            for suffix in 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split():
                pins[spec['index']+'.'+suffix] = '0'*64
        pins.update({str(Path(__file__).resolve()): '0'*64, str(HERE/'analyze16.py'): '0'*64})
        return dict(format='brunello-expanded-benchmark-v1', limits=dict(worker.LIMITS),
                    configs=dict(worker.CONFIGS), guide_count=16, full_run_authorized=False,
                    assemblies=assemblies, pins=pins)

    def test_valid_plan(self):
        with tempfile.TemporaryDirectory() as tmp:
            worker.validate(self.plan(tmp))

    def test_config_scope_and_authorization_drift_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            plan = self.plan(tmp)
            worker.validate(plan)
            for mutation in [dict(worker.CONFIGS, **{'candidate-B': [5, 1, 1]}),
                             dict(worker.CONFIGS, **{'candidate-C': [4, 1, 1]}),
                             dict(worker.CONFIGS, **{'historical': [2, 1, 2]})]:
                with self.assertRaises(ValueError):
                    worker.validate({**plan, 'configs': mutation})
            with self.assertRaises(ValueError):
                worker.validate({**plan, 'guide_count': 1000})
            with self.assertRaises(ValueError):
                worker.validate({**plan, 'limits': dict(worker.LIMITS, rss_kib=1)})
            with self.assertRaises(ValueError):
                worker.validate({**plan, 'full_run_authorized': True})
            with self.assertRaises(ValueError):
                worker.validate({**plan, 'format': 'other-v1'})

    def test_binding_reference_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            plan = self.plan(tmp)
            name = 'HG00735_pat'
            (Path(tmp)/name/'binding.json').write_text(json.dumps(dict(reference_sha256='c'*64)))
            with self.assertRaises(ValueError):
                worker.validate(plan)


class SiteTests(unittest.TestCase):
    def paf(self, guide='g', contig='c', strand='+', start='0', end='20', cigar='20=', cfd='1.0000', cb=None):
        row = [guide, '20', '0', '20', strand, contig, '100', start, end, '20', '20', '255',
               'cg:Z:'+cigar, 'nm:i:0', 'ng:i:0', 'bs:i:0']
        if cfd is not None:
            row.append('cf:f:'+cfd)
        if cb is not None:
            row += ['cb:f:'+cb, 'bu:Z:scored']
        else:
            row.append('bu:Z:test_unavailable')
        return '\t'.join(row)+'\n'

    def test_anchor_key_and_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'x.paf'
            p.write_text(self.paf()+self.paf(strand='-', start='100', end='120'))
            records = analyze.site_records(p)
            self.assertEqual(len(records), 2)
            self.assertIn(('g', 'c', '+', '20'), records)
            self.assertIn(('g', 'c', '-', '100'), records)
            with self.assertRaises(ValueError):
                p.write_text(self.paf()+self.paf())
                analyze.site_records(p)

    def test_nesting_and_representative_changes(self):
        narrow = {('g', 'c', '+', '20'): dict(span=['0', '20'], cigar='20=', cfd='1.0000', cb='0.5', status='scored'),
                  ('g', 'c', '-', '100'): dict(span=['100', '120'], cigar='1X19=', cfd='0.5', cb=None, status='x')}
        broad = dict(narrow)
        broad[('g', 'c', '+', '20')] = dict(span=['0', '20'], cigar='19=1X', cfd='0.9', cb='0.4', status='scored')
        broad[('g', 'c', '+', '50')] = dict(span=['40', '60'], cigar='20=', cfd='0.7', cb='0.2', status='scored')
        result = analyze.nested_comparison({'historical': narrow, 'intermediate': broad,
                                            'candidate-A': broad, 'candidate-B': broad})
        step = result[('historical', 'intermediate')]
        self.assertEqual(step['added_sites'], 1)
        self.assertEqual(step['removed_sites'], 0)
        self.assertEqual(step['shared_representative_differences'], 1)
        self.assertEqual(step['shared_stored_cfd_differences'], 1)
        self.assertEqual(step['shared_neural_availability_changes'], 0)

    def test_cfd_profile_separates_unavailable(self):
        profile = analyze.cfd_profile([None, 0.5, 0.05, 0.004, 0.25])
        self.assertEqual(profile['total'], 5)
        self.assertEqual(profile['defined'], 4)
        self.assertEqual(profile['unavailable'], 1)
        self.assertEqual(profile['ge_0.20'], 2)
        self.assertEqual(profile['ge_0.05'], 3)
        self.assertEqual(profile['max'], 0.5)

    def test_zero_is_not_unavailable(self):
        profile = analyze.cfd_profile([0.0])
        self.assertEqual(profile['defined'], 1)
        self.assertEqual(profile['unavailable'], 0)
        self.assertEqual(profile['ge_0.02'], 0)


if __name__ == '__main__':
    unittest.main()
