import csv
import json
from pathlib import Path
import tempfile
import unittest

import inventory
import collect
import compare_sbds
import qc
import worker


class InventoryTests(unittest.TestCase):
    def test_explicit_assembly_aliases(self):
        self.assertEqual(inventory.historical_id('/x/HG00735_pat_hprc_r2_v1.0.1.fa.gz'), 'HG00735_pat')
        self.assertEqual(inventory.historical_id('/x/hg002v1.1.mat_MT.PanSN.fa.gz'), 'hg002v1.1.mat_MT.PanSN')
        self.assertEqual(inventory.historical_id('/x/chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz'), 'chm13v2.0_maskedY_rCRS.fa.PanSN')

    def test_library_preserves_ids_and_rejects_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'a.fa'
            p.write_text('>original_id\n'+'A'*20+'\n>other_id\n'+'A'*20+'\n')
            self.assertEqual(len(inventory.library(p)), 2)  # identical sequences are not silently deduplicated
            p.write_text('>id\n'+'A'*20+'\n>id\n'+'C'*20+'\n')
            with self.assertRaises(ValueError):
                inventory.library(p)
            p.write_text('>has whitespace\n'+'A'*20+'\n')
            with self.assertRaises(ValueError):
                inventory.library(p)

    def test_strata_are_deterministic_and_input_order_independent(self):
        rows = [('g'+str(n), 'G'*(n % 21)+'A'*(20-n % 21)) for n in range(100)]
        self.assertEqual(inventory.choose_pilot(rows), inventory.choose_pilot(rows[::-1]))


class LegacyReplayTests(unittest.TestCase):
    def test_original_cigar_replay_without_geometry_repair(self):
        result = compare_sbds.replay('A'*20, 'A'*19, '1I19=')
        self.assertEqual(result[3:], (19, 0, 1, 1))
        self.assertEqual(''.join(result[1]), 'A'*20)
        self.assertEqual(''.join(result[2]), '-'+'A'*19)

    def test_bad_base_assertion_is_not_an_eligible_witness(self):
        with self.assertRaises(ValueError):
            compare_sbds.replay('A'*20, 'C'+'A'*19, '20=')
        with self.assertRaises(ValueError):
            compare_sbds.replay('A'*20, 'A'*20, '20M')


class PlanTests(unittest.TestCase):
    def plan(self):
        return dict(format='hprcv2-guarded-pilot-v1', assembly='HG00735_pat', project='sbds',
                    root='/tmp/sbds_pam_aware_v1', full_brunello_authorized=False, fraction=.75, pam='GG', guide_count=1,
                    configurations={'primary': [5, 1, 2], 'historical-domain-control': [2, 1, 1]},
                    reference_sha256='bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b',
                    limits=dict(rss_kib=50331648, min_start_available_kib=67108864,
                        min_running_available_kib=8388608, min_scratch_bytes=60000000000,
                        min_running_disk_bytes=5000000000, stage_timeout_seconds=6600, file_limit_bytes=10737418240),
                    pins={str(Path(worker.__file__).resolve()): 'a'*64})

    def test_valid_safety_plan(self):
        worker.validate(self.plan())

    def test_no_full_library_or_cross_project_authorization(self):
        for key, val in [('full_brunello_authorized', True), ('guide_count', 77441),
                         ('root', '/tmp/brunello_pam_aware_v1'), ('assembly', 'unbound')]:
            plan = self.plan()
            plan[key] = val
            with self.assertRaises(ValueError):
                worker.validate(plan)

    def test_eligibility_and_guard_drift_rejected(self):
        for key, val in [('fraction', .5), ('pam', 'AG'), ('configurations', {'primary': [6, 1, 2]})]:
            plan = self.plan()
            plan[key] = val
            with self.assertRaises(ValueError):
                worker.validate(plan)
        plan = self.plan()
        plan['limits']['rss_kib'] += 1
        with self.assertRaises(ValueError):
            worker.validate(plan)


class CollectionTests(unittest.TestCase):
    def test_pending_is_not_zero_or_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = collect.collect(Path(tmp), 'brunello')
            self.assertEqual(result['status'], 'PENDING_OR_INCOMPLETE')
            self.assertNotIn('unvalidated_linear_sensitivity_scenarios', result)
            self.assertNotIn('eligible_hits', result['configurations']['A'])

    def test_compression_measures_real_stream(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'x'
            p.write_bytes(b'ACGT'*1000)
            self.assertGreater(collect.compressed_bytes(p), 0)
            self.assertLess(collect.compressed_bytes(p), p.stat().st_size)

    def test_search_signature_changes_with_science_not_qc_worker(self):
        p = PlanTests().plan()
        p.update(cli='/cli', native='/native', index='/idx', source='/source', guides='/g')
        p['pins'].update({'/cli': '1'*64, '/native': '2'*64, '/idx.crispr-binding.json': '3'*64,
                          '/g': '4'*64, '/source/native.py': '5'*64})
        original = worker.search_signature(p, 'primary')
        p['pins'][str(Path(worker.__file__).resolve())] = '6'*64
        self.assertEqual(worker.search_signature(p, 'primary'), original)
        p['pins']['/g'] = '7'*64
        self.assertNotEqual(worker.search_signature(p, 'primary'), original)


class AccountingTests(unittest.TestCase):
    def run_fixture(self, cb='0', cfd='1.0000', independent='1.0', mutate=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = [root/n for n in ['raw', 'annotated', 'annotations', 'audit']]
            raw = '\t'.join(['g', '20', '0', '20', '+', 'chr', '23', '0', '20', '20', '20', '255',
                             'cg:Z:20=', 'nm:i:0', 'ng:i:0', 'bs:i:0', 'as:i:0'] + ([] if cfd is None else ['cf:f:'+cfd]))+'\n'
            status = 'scored' if cb is not None else 'test_unavailable_reason'
            annotated = raw.rstrip('\n')+('' if cb is None else '\tcb:f:'+cb)+'\tbu:Z:'+status+'\n'
            ann = dict(row_index='0', guide_id='g', contig='chr', strand='+', start_0based='0', end_0based='20',
                       cigar='20=', pam='AGG', target_sequence='A'*20, aligned_guide='A'*20, aligned_target='A'*20,
                       crispr_bulge_reason=status, crispr_bulge='N/A' if cb is None else cb,
                       cfd='N/A' if cfd is None else cfd, cfd_reason='upstream_cfd_unavailable' if cfd is None else 'available')
            audit = dict(line='1', pam='AGG', target='A'*20, aligned_guide='A'*20, aligned_target='A'*20,
                         mismatches='0', bulges='0', independent_cfd=independent)
            if mutate:
                annotated, ann = mutate(annotated, ann)
            paths[0].write_text(raw)
            paths[1].write_text(annotated)
            for path, row in zip(paths[2:], [ann, audit]):
                with path.open('w') as f:
                    w = csv.DictWriter(f, fieldnames=list(row), delimiter='\t')
                    w.writeheader()
                    w.writerow(row)
            return qc.summarize(*paths, ['g'])

    def test_actual_neural_zero_is_available(self):
        result = self.run_fixture()
        self.assertEqual(result['per_guide']['g']['cb_defined'], 1)
        self.assertEqual(result['per_guide']['g']['view_0.20_rows'], 1)

    def test_unavailable_scores_never_become_zero_or_discarded(self):
        result = self.run_fixture(cb=None, cfd=None, independent='None')
        g = result['per_guide']['g']
        self.assertEqual(g['cb_unavailable'], 1)
        self.assertEqual(g['view_0.20_unscored_preserved'], 1)
        self.assertEqual(g['view_0.20_rows'], 1)

    def test_threshold_uses_full_precision_not_rounded_stored_score(self):
        result = self.run_fixture(cfd='0.0500', independent='0.0499996')
        self.assertEqual(result['per_guide']['g']['view_none_rows'], 1)
        self.assertEqual(result['per_guide']['g']['view_0.05_rows'], 0)

    def test_changed_paf_fields_or_fabricated_unavailable_rejected(self):
        with self.assertRaises(ValueError):
            self.run_fixture(mutate=lambda text, ann: (text.replace('cf:f:1.0000', 'cf:f:0.5'), ann))
        def fabricate(text, ann):
            ann['crispr_bulge'] = '0'
            return text, ann
        with self.assertRaises(ValueError):
            self.run_fixture(cb=None, mutate=fabricate)


if __name__ == '__main__':
    unittest.main()
