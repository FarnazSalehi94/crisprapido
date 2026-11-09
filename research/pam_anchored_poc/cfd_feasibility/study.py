#!/usr/bin/env python3
"""Offline feasibility model, NOT a CRISPRapido backend or FM-index implementation.

Direct PAM enumeration + forward gap templates form the independent oracle.
Reverse traversal groups explicit candidate occurrences instead of SA intervals.
No production code is changed. Run run_study.py for the bounded comparison.
"""
from dataclasses import dataclass
from itertools import groupby, product
import math
from pathlib import Path
import re
import struct
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts/release_validation'))
from audit_emitted import independent_cfd, matrices  # independent existing oracle

DNA = 'ACGT'
RC_TABLE = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(RC_TABLE)[::-1]


def load_tables():
    mm, pm = matrices(ROOT)
    # Production accepts arbitrary user matrices. A future pruning mode must
    # explicitly validate this stronger contract, not silently assume it.
    validate_tables(mm, pm)
    return mm, pm


def validate_tables(mm, pm):
    if not all(math.isfinite(v) and 0 <= v <= 1 for v in (*mm.values(), *pm.values())):
        raise ValueError('pruning requires finite factors in [0,1]')
    required = {f'r{g.replace("T", "U")}:d{rc(t)},{p}'
                for p in range(1, 21) for g in DNA for t in DNA if g != t}
    required |= {f'r{g.replace("T", "U")}:d-,{p}' for p in range(2, 21) for g in DNA}
    required |= {f'r-:d{rc(t)},{p}' for p in range(2, 21) for t in DNA}
    if not required <= mm.keys() or not {''.join(p) for p in product(DNA, repeat=2)} <= pm.keys():
        raise ValueError('pruning requires complete mismatch/gap/PAM factors')


def threshold(value):
    t = float(value)
    if not math.isfinite(t) or not 0 <= t <= 1:
        raise ValueError('research threshold must be finite and in [0,1]')
    return t


def factor(g, t, position, mm):
    if g == t or (position == 1 and '-' in (g, t)):
        return 1.0
    return mm.get(f'r{g.replace("T", "U")}:d{rc(t)},{position}')


def multiply_up(a, b):
    return math.nextafter(a * b, math.inf)


def numeric_upper(product_upper):
    # In addition to directed-up reverse accumulation, cover forward scorer
    # rounding (<=21 products), with ample relative and subnormal allowance.
    return math.nextafter(product_upper * (1 + 64 * sys.float_info.epsilon)
                          + 64 * math.ulp(0.0), math.inf)


def f32(x):
    return struct.unpack('f', struct.pack('f', x))[0]


def passes(matches, length, frac):
    return matches > 0 and f32(f32(matches / length) * 100) >= f32(f32(frac) * 100)


@dataclass(frozen=True)
class Alignment:
    key: tuple  # contig, strand, PAM-adjacent target boundary
    start: int  # genomic target start, not oriented offset
    end: int
    target: str
    cigar: str
    ops: str
    mm: int
    ins: int
    dele: int
    score: object

    def order(self):
        return self.mm + self.ins + self.dele, self.mm, self.cigar, self.start


def alignment(guide, target, ops, key, start, end, pam, tables, m, frac):
    qi = ti = mismatches = matches = 0
    aq, at, true_ops = [], [], []
    for op in ops:
        if op == 'I':
            g, t = guide[qi], '-'
            qi += 1
        elif op == 'D':
            g, t = '-', target[ti]
            ti += 1
        else:
            g, t = guide[qi], target[ti]
            qi += 1
            ti += 1
            matches += g == t
            mismatches += g != t
            op = '=' if g == t else 'X'
        aq.append(g)
        at.append(t)
        true_ops.append(op)
    assert qi == len(guide) and ti == len(target)
    if mismatches > m or not passes(matches, len(guide), frac):
        return None
    grouped = [(len(list(xs)), op) for op, xs in groupby(true_ops)]
    cigar = ''.join(f'{n}{op}' for n, op in grouped)
    score = independent_cfd(''.join(aq), ''.join(at), grouped, pam, *tables)
    return Alignment(key, start, end, target, cigar, ''.join(true_ops), mismatches,
                     ops.count('I'), ops.count('D'), score)


def templates(length, z):
    """Forward all one-event gap words; independent of reverse search state."""
    yield '=' * length
    for size in range(1, z + 1):
        for word in map(''.join, product('ID', repeat=size)):
            for prefix in range(length - word.count('I') + 1):
                yield '=' * prefix + word + '=' * (length - prefix - word.count('I'))


def occurrences(refs, pam):
    out = []
    for contig, seq in refs.items():
        for strand in '+-':
            oriented = seq if strand == '+' else rc(seq)
            for end in range(1, len(oriented) - 2):
                if oriented[end] in DNA and oriented[end+1:end+3] == pam:
                    anchor = end if strand == '+' else len(seq) - end
                    out.append((contig, strand, anchor, oriented, end))
    return out


def oracle(refs, guide, pam, m, z, frac, tables):
    out = set()
    for contig, strand, anchor, oriented, end in occurrences(refs, pam):
        for ops in templates(len(guide), z):
            span = len(ops) - ops.count('I')
            if span <= 0 or span > end:
                continue
            target = oriented[end-span:end]
            lo, hi = ((end-span, end) if strand == '+' else
                      (len(oriented)-end, len(oriented)-end+span))
            a = alignment(guide, target, ops, (contig, strand, anchor), lo, hi,
                          oriented[end:end+3], tables, m, frac)
            if a:
                out.add(a)
    return out


def representatives(alignments):
    result = {}
    for a in alignments:
        if a.key not in result or a.order() < result[a.key].order():
            result[a.key] = a
    return result


def keep(a, tau):
    return a.score is None or a.score >= tau


def may_be_unavailable(length, has_n, i, state, gaplen, ins, dele, z):
    # Conservative global N fallback; the real index does not annotate N
    # positions in SA intervals. Mixed I/D has 20+D columns and no projection.
    if length != 20 or has_n or (ins and dele):
        return True
    if state == 2 or z < 2:
        return False
    if state == 0:
        return i >= 0  # an I and D can still be introduced
    if gaplen >= z:
        return False
    return bool(ins or (dele and i >= 0))


def reverse_model(refs, guide, pam, m, z, frac, tables, tau=None):
    """Explicit occurrence-set search; state count is NOT a Columba measurement.

    N->A shadow models an admissible substituted index (not builder randomness).
    Emission always replays the true reference. Global N presence disables all
    score pruning: conservative and deliberately not an N-aware index feature.
    """
    mmatrix, pmatrix = tables
    seeds = occurrences({c: s.replace('N', 'A') for c, s in refs.items()}, pam)
    has_n = any('N' in s for s in refs.values())
    stats = dict(states=0, pruned=0, emitted_paths=0)
    out = set()

    def walk(candidates, i, mm, state, gaplen, ins, dele, span, matches, revops, upper, columns):
        stats['states'] += 1
        if mm > m or not passes(matches + i + 1, len(guide), frac):
            return
        if (tau is not None and not may_be_unavailable(len(guide), has_n, i, state, gaplen, ins, dele, z)
                and numeric_upper(upper) < tau):
            stats['pruned'] += 1
            return
        if i < 0:
            for c, strand, anchor, _, end in candidates:
                if span <= 0:
                    continue
                true_ori = refs[c] if strand == '+' else rc(refs[c])
                observed_pam = true_ori[end:end+3]
                if observed_pam[0] not in DNA or observed_pam[1:] != pam:
                    continue
                lo, hi = ((end-span, end) if strand == '+' else
                          (len(true_ori)-end, len(true_ori)-end+span))
                a = alignment(guide, true_ori[end-span:end], revops[::-1], (c, strand, anchor),
                              lo, hi, observed_pam, tables, m, frac)
                if a:
                    stats['emitted_paths'] += 1
                    out.add(a)
        # PAM-proximal columns have invariant scored positions 20,19,...,1.
        # Pure D projection retains the last 20 columns; distal excess is ignored.
        def updated(g, t):
            f = factor(g, t, 20-columns, mmatrix) if columns < 20 else 1.0
            return multiply_up(upper, f) if f is not None else upper

        groups = {c: tuple(x for x in candidates if x[4] > span and x[3][x[4]-span-1] == c)
                  for c in DNA}
        if state != 2 and gaplen < z:
            for c, subset in groups.items():
                if subset:
                    walk(subset, i, mm, 1, gaplen+1, ins, dele+1, span+1, matches,
                         revops+'D', updated('-', c), columns+1)
        if i < 0:
            return
        if state != 2 and gaplen < z:
            walk(candidates, i-1, mm, 1, gaplen+1, ins+1, dele, span, matches,
                 revops+'I', updated(guide[i], '-'), columns+1)
        for c, subset in groups.items():
            if subset:
                mismatch = c != guide[i]
                walk(subset, i-1, mm+mismatch, 2 if state else 0, gaplen, ins, dele,
                     span+1, matches+(not mismatch), revops+('X' if mismatch else '='),
                     updated(guide[i], c), columns+1)

    for n in DNA:
        subset = tuple(x for x in seeds if x[3][x[4]] == n)
        if subset:
            walk(subset, len(guide)-1, 0, 0, 0, 0, 0, 0, 0, '', pmatrix[pam], 0)
    return out, stats


def fixtures():
    """Finite exhaustive suffix variations + explicit gap/N controls, both strands."""
    guide = 'GAAACAGTCGATTTTATCAC'
    refs = {}
    def add(name, target, pam='GG'):
        for strand in '+-':
            s = target + 'A' + pam
            refs[name + strand] = s if strand == '+' else rc(s)
    for suffix in map(''.join, product(DNA, repeat=4)):
        add('suffix_' + suffix, guide[:-4] + suffix)
    for k in (0, 1, 5, 14, 19, 20):
        for n in (1, 2):
            add(f'D{k}_{n}', guide[:k] + 'T'*n + guide[k:])
            if k+n <= 20:
                add(f'I{k}_{n}', guide[:k] + guide[k+n:])
        if k < 20:
            add(f'mixed{k}', guide[:k] + ('C' if guide[k] != 'C' else 'A') + guide[k+1:])
    add('repeat', 'AC'*11)
    # Non-N references for all 16 PAM suffixes, also zero-PAM-score controls.
    for pam in map(''.join, product(DNA, repeat=2)):
        add('pam_' + pam, guide, pam)
    nrefs = dict(refs)
    add_n = [('targetN', guide[:16]+'N'+guide[17:]+'AGG'),
             ('distalN_D', 'NN'+guide+'AGG'),
             ('pamN', guide+'NGG'), ('suffixN', guide+'AGN')]
    for name, s in add_n:
        nrefs[name+'+'] = s
        nrefs[name+'-'] = rc(s)
    return guide, refs, nrefs
