#!/usr/bin/env python3
"""Experimental backend runner: search a private, verified reference/index snapshot."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'compatibility'))
from index_manifest import SUFFIXES, verify


def main():
    p = argparse.ArgumentParser()
    for name in ('binary', 'reference', 'index', 'mismatches', 'bulges', 'size', 'fraction', 'pam'):
        p.add_argument('--' + name, required=True)
    a = p.parse_args()
    with tempfile.TemporaryDirectory(prefix='crisprapido-pam-') as td:
        root = Path(td)
        ref, idx, guides = root / 'reference.fa', root / 'index', root / 'guides.fa'
        shutil.copyfile(a.reference, ref)
        for suffix in (*SUFFIXES, 'crispr-binding.json'):
            shutil.copyfile(a.index + '.' + suffix, str(idx) + '.' + suffix)
        guides.write_text(sys.stdin.read(), encoding='ascii')
        verify(ref, idx)
        for path in root.iterdir():
            path.chmod(0o400)
        run = subprocess.run([a.binary, str(idx), str(ref), str(guides), a.mismatches,
                              a.size, a.fraction, a.pam, a.bulges], capture_output=True, text=True)
        sys.stderr.write(run.stderr)
        if run.returncode:
            raise RuntimeError('PAM search executable failed with status ' + str(run.returncode))
        verify(ref, idx)
        # Extract targets only from the exact verified snapshot used by the search.
        seqs = {}
        name = None
        for line in ref.read_text().splitlines():
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in seqs:
                    raise ValueError('duplicate reference record')
                seqs[name] = ''
            elif name is not None:
                seqs[name] += ''.join(line.split()).upper()
        lines = run.stdout.splitlines()
        header = 'guide\tcontig\tstrand\tstart\tend\tanchor\tpam\tcigar\tmismatches\tinserted_guide_bases\tdeleted_guide_bases\tbulge_events\tmatch_count\tmatch_fraction'
        try:
            start = lines.index(header)
        except ValueError:
            raise ValueError('search output header missing')
        print(header + '\ttarget\tcontig_length')
        complement = str.maketrans('ACGTN', 'TGCAN')
        for line in lines[start + 1:]:
            fields = line.split('\t')
            if len(fields) != 14:
                raise ValueError('invalid search output')
            seq = seqs[fields[1]]
            lo, hi = int(fields[3]), int(fields[4])
            if not 0 <= lo <= hi <= len(seq):
                raise ValueError('invalid search coordinates')
            target = seq[lo:hi]
            if fields[2] == '-':
                target = target.translate(complement)[::-1]
            elif fields[2] != '+':
                raise ValueError('invalid search strand')
            print(line + '\t' + target + '\t' + str(len(seq)))

if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, RuntimeError) as e:
        print('experimental PAM backend: ' + str(e), file=sys.stderr)
        sys.exit(1)
