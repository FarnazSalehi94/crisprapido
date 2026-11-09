#!/usr/bin/env python3
"""Regression tests for offline bounds/semantics, not production search tests."""
from fractions import Fraction
import math
import unittest

from study import (alignment, factor, fixtures, keep, load_tables, may_be_unavailable,
                   multiply_up, numeric_upper, oracle, representatives, reverse_model,
                   rc, templates, threshold, validate_tables)


class FeasibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables = load_tables()
        cls.guide = 'GAAACAGTCGATTTTATCAC'

    def test_matrix_contract(self):
        mm, pm = self.tables
        for value in (math.nan, math.inf, -0.1, 1.01):
            bad = dict(mm)
            bad[next(iter(bad))] = value
            with self.assertRaises(ValueError):
                validate_tables(bad, pm)
        bad = dict(mm)
        del bad['rA:dG,20']
        with self.assertRaises(ValueError):
            validate_tables(bad, pm)

    def test_threshold_validation(self):
        for x in (-1, 1.1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                threshold(x)
        self.assertEqual(threshold(0), 0)
        self.assertEqual(threshold(1), 1)

    def test_supported_and_unsupported_bulges(self):
        g = self.guide
        cases = [(g, '='*20, True), (g[:5]+g[6:], '='*5+'I'+'='*14, True),
                 ('T'+g, 'D'+'='*20, True), ('TT'+g, 'DD'+'='*20, True),
                 (g, 'DI'+'='*19, False), (g, 'ID'+'='*19, False),
                 ('N'+g[1:], 'X'+'='*19, False)]
        for t, ops, supported in cases:
            a = alignment(g, t, ops, ('c','+',len(t)), 0,len(t),'AGG',self.tables,20,0)
            self.assertEqual(a.score is not None, supported)
            if not supported:
                self.assertTrue(keep(a, 1))
        # A reference N can be projected away and leave a defined score.
        a = alignment(g, 'N'+g, 'D'+'='*20, ('c','+',21), 0,21,'AGG',self.tables,20,0)
        self.assertEqual(a.score, 1)

    def test_projected_suffix_bound(self):
        g = self.guide
        targets = [g, 'T'+g, g[:5]+'TT'+g[5:], g[:-1], 'T'+g[1:]]
        checked = 0
        for ops in templates(20, 2):
            span = len(ops)-ops.count('I')
            for t in targets:
                if len(t) != span:
                    continue
                a = alignment(g,t,ops,('c','+',span),0,span,'AGG',self.tables,20,0)
                if a is None or a.score is None:
                    continue
                qi = ti = 0
                cols = []
                for op in ops:
                    if op == 'I':
                        cols.append((g[qi], '-')); qi += 1
                    elif op == 'D':
                        cols.append(('-', t[ti])); ti += 1
                    else:
                        cols.append((g[qi],t[ti])); qi += 1; ti += 1
                u = 1.0
                exact = Fraction(1)
                for j, (gb,tb) in enumerate(reversed(cols), 1):
                    f = factor(gb,tb,21-j,self.tables[0]) if j<=20 else 1.0
                    self.assertIsNotNone(f)
                    exact *= Fraction(f)
                    u = multiply_up(u,f)
                    self.assertGreaterEqual(Fraction(u), exact)
                    self.assertGreaterEqual(numeric_upper(u), a.score)
                    checked += 1
        self.assertGreater(checked, 1000)

    def test_unknown_guards(self):
        self.assertTrue(may_be_unavailable(20,False,19,0,0,0,0,2))
        self.assertTrue(may_be_unavailable(20,False,18,1,1,1,0,2))
        self.assertTrue(may_be_unavailable(20,False,18,2,2,1,1,2))
        self.assertFalse(may_be_unavailable(20,False,18,2,1,1,0,2))
        self.assertTrue(may_be_unavailable(20,True,18,2,1,1,0,2))
        self.assertTrue(may_be_unavailable(19,False,18,0,0,0,0,0))

    def test_bounded_oracle_and_exact_threshold(self):
        g = self.guide
        refs = {'exact':g+'AGG', 'mismatch':g[:-1]+'AAGG', 'mixed':g[:5]+'C'+g[6:]+'AGG',
                'insertion':g[:7]+g[8:]+'AGG', 'deletion':g[:8]+'T'+g[8:]+'AGG'}
        for m,z,f in [(0,0,.75),(5,1,.75),(8,2,.5),(20,2,0)]:
            truth = oracle(refs,g,'GG',m,z,f,self.tables)
            base, _ = reverse_model(refs,g,'GG',m,z,f,self.tables)
            self.assertEqual(base,truth)
            # At equality a defined score must be retained, not rounded away.
            thresholds = {0,1,.1} | {x.score for x in truth if x.score is not None}
            for tau in thresholds:
                found, _ = reverse_model(refs,g,'GG',m,z,f,self.tables,tau)
                self.assertLessEqual({x for x in truth if keep(x,tau)},found)
                keys = {x.key for x in found}
                final = {k:v for k,v in representatives(x for x in truth if x.key in keys).items() if keep(v,tau)}
                expected = {k:v for k,v in representatives(truth).items() if keep(v,tau)}
                self.assertEqual(final,expected)

    def test_representative_refinement_needed(self):
        g = self.guide
        refs = {'I14_1-': rc(g[:14]+g[15:]+'AGG')}
        truth = oracle(refs,g,'GG',5,1,.75,self.tables)
        canonical = representatives(truth)
        discovered, _ = reverse_model(refs,g,'GG',5,1,.75,self.tables,.1)
        expected = {k:v for k,v in canonical.items() if keep(v,.1)}
        naive = {k:v for k,v in representatives(discovered).items() if keep(v,.1)}
        self.assertNotEqual(naive,expected)
        keys = {x.key for x in discovered}
        refined = {k:v for k,v in representatives(x for x in truth if x.key in keys).items() if keep(v,.1)}
        self.assertEqual(refined,expected)

    def test_unavailable_mixed_representative_survives(self):
        g = self.guide
        refs = {'mixed_only': g[:5]+'C'+g[6:]+'AGG'}
        truth = oracle(refs,g,'GG',0,2,.75,self.tables)
        canonical = representatives(truth)
        self.assertTrue(canonical)
        self.assertTrue(all(a.score is None for a in canonical.values()))
        found, _ = reverse_model(refs,g,'GG',0,2,.75,self.tables,1)
        self.assertLessEqual(truth,found)
        self.assertEqual(representatives(found),canonical)

    def test_n_fallback_and_non20_guides(self):
        for g in (self.guide, self.guide[:19]):
            refs = {'n':'N'+g[1:]+'AGG','exact':g+'AGG'}
            baseline, before = reverse_model(refs,g,'GG',20,2,0,self.tables)
            pruned, after = reverse_model(refs,g,'GG',20,2,0,self.tables,1)
            self.assertEqual(baseline,pruned)
            self.assertEqual(before,after)


if __name__ == '__main__':
    unittest.main()
