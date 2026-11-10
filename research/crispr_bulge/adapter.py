#!/usr/bin/env python3
"""Stdlib-only, alignment-preserving input conversion. No search or rescoring CFD."""
import hashlib
from pathlib import Path
import re

DNA = set('ACGT')
IUPAC = set('ACGTRYSWKMBDHVN')
COMP = str.maketrans('ACGTRYMKSWBDHVN', 'TGCAYRKMSWVHDBN')


def rc(s):
    return s.translate(COMP)[::-1]


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


class IndexedFasta:
    """Build an in-memory FAI dictionary, not an in-memory reference sequence.

    Plain, regularly wrapped FASTA only. The last line may be short. Scan hashes
    the exact bytes and validates alphabet/names/wrapping. No input file is modified.
    """
    def __init__(self, path):
        self.path = Path(path)
        self.rows = {}
        sha = hashlib.sha256()
        name = None
        length = offset = bases = width = 0
        short = False
        with self.path.open('rb') as f:
            while True:
                position = f.tell()
                line = f.readline()
                if not line:
                    break
                sha.update(line)
                if line.startswith(b'>'):
                    if name is not None:
                        if not length:
                            raise ValueError('empty reference contig')
                        self.rows[name] = (length, offset, bases, width)
                    fields = line[1:].split()
                    if not fields:
                        raise ValueError('empty reference ID')
                    name = fields[0].decode('ascii')
                    if name in self.rows:
                        raise ValueError('duplicate reference ID')
                    length = bases = width = 0
                    short = False
                else:
                    seq = line.rstrip(b'\r\n')
                    if name is None or not seq or not set(seq.decode('ascii').upper()) <= IUPAC:
                        raise ValueError('unsupported FASTA syntax/alphabet')
                    if short:
                        raise ValueError('irregular FASTA wrapping')
                    if not bases:
                        offset, bases, width = position, len(seq), len(line)
                    elif len(seq) > bases or (len(seq) == bases and len(line) != width):
                        # A final full-width line without its newline is valid.
                        if len(seq) == bases and not line.endswith(b'\n'):
                            short = True
                        else:
                            raise ValueError('irregular FASTA wrapping')
                    elif len(seq) < bases:
                        short = True
                    length += len(seq)
        if name is None or not length:
            raise ValueError('empty reference')
        self.rows[name] = (length, offset, bases, width)
        self.sha256 = sha.hexdigest()
        self.file = self.path.open('rb')

    def fetch(self, name, lo, hi):
        length, offset, bases, width = self.rows[name]
        if not 0 <= lo <= hi <= length:
            raise ValueError('reference interval out of bounds')
        pieces = []
        while lo < hi:
            n = min(hi-lo, bases-lo % bases)
            self.file.seek(offset + (lo // bases)*width + lo % bases)
            part = self.file.read(n).decode('ascii').upper()
            if len(part) != n or not set(part) <= IUPAC:
                raise ValueError('reference changed or invalid FASTA offsets')
            pieces.append(part)
            lo += n
        return ''.join(pieces)

    def close(self):
        self.file.close()


def guides_from_fasta(path):
    guides = {}
    name = None
    for line in Path(path).read_text(encoding='ascii').splitlines():
        if line.startswith('>'):
            parts = line[1:].split()
            if not parts or parts[0] in guides:
                raise ValueError('empty/duplicate guide ID')
            name = parts[0]
            guides[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError('sequence before guide header')
            guides[name] += ''.join(line.split()).upper()
    if not guides or any(not s for s in guides.values()):
        raise ValueError('empty guides')
    return guides


def parse_paf(raw):
    fields = raw.rstrip('\r\n').split('\t')
    if len(fields) < 12 or fields[4] not in ('+', '-'):
        raise ValueError('invalid PAF row')
    tags = {}
    for text in fields[12:]:
        pieces = text.split(':', 2)
        if len(pieces) != 3 or len(pieces[0]) != 2 or pieces[0] in tags:
            raise ValueError('invalid/duplicate PAF tag')
        tags[pieces[0]] = (pieces[1], pieces[2])
    if 'cb' in tags or 'bu' in tags:
        raise ValueError('input already contains CRISPR-Bulge annotations')
    return fields, tags


def convert(raw, guides, reference):
    f, tags = parse_paf(raw)
    g = guides[f[0]]
    qlen, qlo, qhi = map(int, f[1:4])
    lo, hi = int(f[7]), int(f[8])
    if qlen != len(g) or int(f[6]) != reference.rows[f[5]][0]:
        raise ValueError('guide/reference length mismatch')
    if not 0 <= qlo <= qhi <= qlen or not 0 <= lo <= hi <= int(f[6]):
        raise ValueError('invalid PAF interval')
    key = f[0], f[5], f[4], hi if f[4] == '+' else lo
    if qlo != 0 or qhi != qlen:
        return key, 'unsupported_partial_guide', None
    if tags.get('cg', ('', ''))[0] != 'Z':
        raise ValueError('missing CIGAR')
    cigar = tags['cg'][1]
    if not re.fullmatch(r'(?:[1-9][0-9]*[=XMID])+', cigar):
        raise ValueError('unsupported/malformed CIGAR')
    ops = [(int(n), op) for n, op in re.findall(r'(\d+)([=XMID])', cigar)]
    if (sum(n for n, op in ops if op != 'D') != len(g) or
            sum(n for n, op in ops if op != 'I') != hi-lo):
        raise ValueError('CIGAR consumption mismatch')
    target = reference.fetch(f[5], lo, hi)
    if f[4] == '-':
        target = rc(target)
    qi = ti = mismatch = ins = dele = events = 0
    previous_gap = False
    ag, at = [], []
    for n, op in ops:
        gap = op in 'ID'
        events += gap and not previous_gap
        previous_gap = gap
        for _ in range(n):
            if op == 'I':
                a, b = g[qi], '-'
                qi += 1
                ins += 1
            elif op == 'D':
                a, b = '-', target[ti]
                ti += 1
                dele += 1
            else:
                a, b = g[qi], target[ti]
                qi += 1
                ti += 1
                if a in DNA and b in DNA:
                    if (op == '=' and a != b) or (op == 'X' and a == b):
                        raise ValueError('CIGAR/reference assertion failed')
                mismatch += a != b
            ag.append(a)
            at.append(b)
    if len(g) != 20:
        return key, 'unsupported_guide_length', None
    if not set(g) <= DNA:
        return key, 'ambiguous_guide', None
    if not set(target) <= DNA:
        return key, 'ambiguous_reference', None
    if (hi+3 > int(f[6]) if f[4] == '+' else lo < 3):
        return key, 'incomplete_pam', None
    pam = reference.fetch(f[5], hi, hi+3) if f[4] == '+' else rc(reference.fetch(f[5], lo-3, lo))
    if not set(pam) <= DNA:
        return key, 'ambiguous_pam', None
    if pam[1:] != 'GG':
        return key, 'unsupported_pam', None
    if len(ag)+3 > 24:
        return key, 'alignment_exceeds_24_columns', None
    if ins and dele:
        return key, 'unsupported_mixed_bulge', None
    if events > 1:
        return key, 'unsupported_multiple_bulge_events', None
    if ins+dele > 1:
        return key, 'unsupported_multibase_bulge', None
    # Conservative documented GUIDE-seq preprocessing domain, not a hard neural
    # network input limit. Never impose these limits on either search backend.
    if mismatch > (4 if ins+dele else 6):
        return key, 'outside_documented_mismatch_domain', None
    return key, 'scored', {'guide': ''.join(ag)+'NGG', 'target': ''.join(at)+pam}


def encode_tokens(rows):
    """Independent categorical encoding for the pinned C_2 checkpoint.

    Alphabet order A,C,G,T,-; left GAP padding, then pair index 5*g+t.
    Structural N in the guide PAM maps to the observed target base, as upstream.
    True ambiguous sequence input must have been rejected before this function.
    """
    out = []
    alphabet = 'ACGT-'
    for row in rows:
        g, t = row['guide'], row['target']
        if len(g) != len(t) or not 23 <= len(g) <= 24 or g[-3:] != 'NGG':
            raise ValueError('invalid model input geometry/PAM')
        if set(g[:-3]) - set(alphabet) or set(t) - set(alphabet) or t[-2:] != 'GG':
            raise ValueError('ambiguous/unsupported model input')
        gap_bases = g[:-3].count('-') + t[:-3].count('-')
        mismatches = sum(a != b for a, b in zip(g[:-3], t[:-3]) if a != '-' and b != '-')
        if (len(g[:-3].replace('-', '')) != 20 or gap_bases > 1 or
                mismatches > (4 if gap_bases else 6) or any(a == b == '-' for a, b in zip(g, t))):
            raise ValueError('outside conservative model domain')
        g, t = g.rjust(24, '-'), t.rjust(24, '-')
        tokens = []
        for a, b in zip(g, t):
            if a == 'N':
                a = b
            tokens.append(5*alphabet.index(a)+alphabet.index(b))
        out.append(tokens)
    return out
