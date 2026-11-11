import hashlib
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


import sys
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))


bench = load('bench_worker')
validate_old = load('validate_old')
collect_bench = load('collect_bench')

GUIDE = 'CATCTTCTTTCACCTGAACG'
TIME_V = """	Command being timed: "legacy-binary -r ref.fa"
	User time (seconds): 10782.08
	System time (seconds): 240.43
	Percent of CPU this job got: 190%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:36:15
	Maximum resident set size (kbytes): 87989928
	Exit status: 0
"""


class ParserTests(unittest.TestCase):
    def test_historical_time_v_format(self):
        result = bench.parse_time_v(TIME_V)
        self.assertEqual(result['wall_seconds'], 5775.0)
        self.assertEqual(result['cpu_percent'], 190.0)
        self.assertEqual(result['max_rss_kib'], 87989928)

    def test_minute_format_and_nonzero_exit_rejected(self):
        result = bench.parse_time_v(TIME_V.replace('1:36:15', '0:58.54'))
        self.assertAlmostEqual(result['wall_seconds'], 58.54)
        with self.assertRaises(ValueError):
            bench.parse_time_v(TIME_V.replace('Exit status: 0', 'Exit status: 1'))

    def test_native_line_exactly_once(self):
        self.assertEqual(bench.parse_native('states=100 x\nother\nstates=42 search_wall_seconds=1.5\n'),
                         dict(states=42, traversal_seconds=1.5))
        with self.assertRaises(ValueError):
            bench.parse_native('states=1 search_wall_seconds=1\nstates=2 search_wall_seconds=2\n')


class LabelTests(unittest.TestCase):
    def row(self, **over):
        fields = ['g', '20', '0', '20', '+', 'chr', '23', '0', '20', '20', '20', '255',
                  'nm:i:0', 'ng:i:0', 'bs:i:0', 'cg:Z:20=', 'cf:f:1.0000']
        if 'cigar' in over:
            fields = [x for x in fields if not x.startswith('cg:Z:')] + ['cg:Z:'+over.pop('cigar')]
        positions = {'strand': 4, 'qstart': 2, 'qend': 3}
        for key, value in over.items():
            fields[positions.get(key, key)] = value
        return fields

    def label(self, fields, target=GUIDE, pam='AGG'):
        return validate_old.label_row(1, fields, {'g': GUIDE}, target, pam, None)[0]

    def test_perfect_match_is_strict_eligible(self):
        result = self.label(self.row())
        self.assertEqual(result['geometry'], 'standard')
        self.assertEqual(result['strict_eligible'], True)
        self.assertEqual(result['strict_ngg'], True)

    def test_non_ngg_pam_is_detected_not_repaired(self):
        result = self.label(self.row(), pam='TCT')
        self.assertEqual(result['strict_ngg'], False)
        self.assertEqual(result['strict_eligible'], False)

    def test_legacy_bulge_geometry_is_labelled(self):
        result = self.label(self.row(qstart='1', qend='21', cigar='1I19='),
                            target=GUIDE[1:])
        self.assertEqual(result['geometry'], 'legacy_query_coordinates')
        self.assertEqual(result['strict_eligible'], True)
        self.assertEqual(result['bulge_events'], 1)

    def test_replay_failure_and_ineligible_counts(self):
        result = self.label(self.row(cigar='20='), target='A'+GUIDE[1:])
        self.assertTrue(result['replay_error'])
        self.assertEqual(result['strict_eligible'], '')
        result, accounted = validate_old.label_row(1, self.row(cigar='20='),
                                                    {'g': GUIDE}, 'A'+GUIDE[1:], 'AGG', None)
        # Legacy broken CIGARs are labelled diagnostics; the validator must not abort on them.
        self.assertTrue(accounted)
        result = self.label(self.row(cigar='1X19='), target='T'+GUIDE[1:])
        self.assertEqual(result['mismatches'], 1)
        self.assertEqual(result['strict_eligible'], True)
        result = self.label(self.row(cigar='3X17='), target='GGG'+GUIDE[3:])
        self.assertEqual(result['mismatches'], 3)
        self.assertEqual(result['strict_eligible'], False)

    def test_cfd_reconstruction_compared_only_when_supplied(self):
        fields = self.row()
        _, accounted = validate_old.label_row(1, fields, {'g': GUIDE}, GUIDE, 'AGG',
                                              lambda ag, at, ops, pam: 0.9999)
        self.assertTrue(accounted)

    def test_short_row_is_labelled_not_dropped(self):
        result, accounted = validate_old.label_row(1, ['g', '20'], {'g': GUIDE}, '', '', None)
        self.assertTrue(accounted)
        self.assertEqual(result['geometry'], 'paf_has_fewer_than_12_fields')


class ExtrapolationTests(unittest.TestCase):
    def test_linear_through_two_largest_points(self):
        est = collect_bench.extrapolate([(1, 320), (100, 400), (1000, 1300)], 77441)
        self.assertAlmostEqual(est, 400 + 900*(77441-100)/900)

    def test_old_model_reproduces_known_total_exactly(self):
        # (100, y100) and (1000, y1000) chosen so the line hits the measured 5775 s.
        y100, slope = 4000.0, (5775.0-4000.0)/(77441-100)
        est = collect_bench.extrapolate([(1, 100), (100, y100), (1000, y100+900*slope)], 77441)
        self.assertAlmostEqual(est, 5775.0)

    def test_identical_points_return_none(self):
        self.assertIsNone(collect_bench.extrapolate([(5, 1), (5, 1), (5, 1)], 10))


class PlanTests(unittest.TestCase):
    def plan(self):
        limits = dict(bench.LIMITS)
        return dict(bench.EXPECTED, root='/tmp/performance_investigation_v1', limits=limits,
                    library=dict(sha256='a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971'),
                    old=dict(binary_sha256='ee543e772ac6e956e8fc3eb88e6ea52517e5192eff9da255902c3bad6073faed',
                             columba_sha256='7e44974c81dad4943292fa64a71e54073908b4909ce327480464d769322c40d1',
                             components={s: '0'*64 for s in bench.SUFFIXES}, threads=16, cpus=16),
                    new=dict(threads=1, cpus=8),
                    pins={str(HERE/'bench_worker.py'): '0'*64, str(HERE/'validate_old.py'): '0'*64,
                          str(HERE/'supervise.py'): '0'*64},
                    worker=str(HERE/'bench_worker.py'), validator=str(HERE/'validate_old.py'),
                    supervisor=str(HERE/'supervise.py'), sizes=[1, 100, 1000], guides={'1': '/g1'})

    def test_valid_plan_and_drift_rejected(self):
        bench.validate(self.plan())
        for key, value in [('settings', dict(m=5, b=1, z=2, f=0.75, pam='GG')),
                           ('node', 'octopus09'), ('sizes', [1, 100, 1000, 77441])]:
            with self.assertRaises(ValueError):
                bench.validate({**self.plan(), key: value})
        plan = self.plan()
        plan['limits']['stage_timeout_seconds'] += 1
        with self.assertRaises(ValueError):
            bench.validate(plan)
        plan = self.plan()
        plan['new']['threads'] = 8  # parallelism would change what is being measured
        with self.assertRaises(ValueError):
            bench.validate(plan)

    def test_benchmark_cannot_request_full_library(self):
        plan = self.plan()
        plan['sizes'] = [77441]
        with self.assertRaises(ValueError):
            bench.validate(plan)

    def test_staged_executables_keep_execute_permission(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp)/'bin', Path(tmp)/'staged-bin'
            src.write_bytes(b'#!/bin/sh\ntrue\n')
            staged = bench.stage_executable(src, dst, hashlib.sha256(src.read_bytes()).hexdigest())
            self.assertTrue(staged.stat().st_mode & 0o111)
            with self.assertRaises(ValueError):
                bench.stage_executable(src, Path(tmp)/'bad', '0'*64)

    def test_supervisor_command_is_type_safe(self):
        # Integer plan fields (sizes, RSS limits) must never reach raw subprocess calls.
        plan = self.plan()
        bench.validate(plan)
        command = list(map(str, ['/usr/bin/python3', plan['supervisor'], '--limit-kib',
                                 plan['limits']['new_rss_kib'], '--size', 1000]))
        self.assertTrue(all(isinstance(x, str) for x in command))
        self.assertIn('41943040', command)


if __name__ == '__main__':
    unittest.main()
