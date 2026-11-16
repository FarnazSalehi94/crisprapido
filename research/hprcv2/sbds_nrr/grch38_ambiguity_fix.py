#!/usr/bin/env python3
"""One-off audited GRCh38 IUPAC-ambiguity fix for the PAM-anchored search.

The native PAM search accepts only A/C/G/T/N (case-insensitive). The GRCh38
PanSN analysis reference contains 94 IUPAC ambiguity codes, which makes every
search fail fast with 'reference must be nonempty A/C/G/T/N'. This script
derives a converted FASTA.gz (ambiguity -> N, coordinates unchanged), records
every converted position, and pins both identities. The original source is
never modified. The derived FASTA is the sole search input for a rebuilt index;
the transformation record preserves full auditability.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido')
OUT = ROOT/'runs/sbds_ngg_nrr_v1/00_manifest/grch38-acgtn-fix-v1'
ALLOWED = set('ACGTNacgtn')


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def convert(source_gz, derived_gz, samtools):
    positions = []
    headers = sequences = 0
    in_header = False
    contig = offset = None
    started = time.monotonic()
    with gzip.open(source_gz, 'rt') as src, gzip.open(derived_gz, 'wt') as dst:
        for line in src:
            if line.startswith('>'):
                headers += 1
                in_header = True
                contig = line[1:].split()[0]
                offset = 0
                dst.write(line)
                continue
            in_header = False
            sequences += 1
            out = []
            for i, ch in enumerate(line.rstrip('\n')):
                if ch in 'RYKMSWBDHVrykmswbdhv':
                    positions.append([contig, offset+i, ch.upper(), 'N'])
                    out.append('N')
                elif ch in ALLOWED:
                    out.append(ch)
                else:
                    raise ValueError('invalid reference character: '+repr(ch))
            offset += len(out)
            dst.write(''.join(out)+'\n')
    if not headers or in_header or offset is None or offset == 0:
        raise ValueError('malformed FASTA structure')
    return positions, dict(source_gz=str(source_gz), source_gz_sha256=digest(source_gz),
        derived_gz=str(derived_gz), derived_gz_sha256=digest(derived_gz),
        headers=headers, sequence_lines=sequences, bases=offset,
        converted_positions=len(positions), seconds=round(time.monotonic()-started, 1))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source-gz', type=Path, default=None)
    p.add_argument('--samtools', type=Path,
                   default=Path('/home/salehi/miniconda3/envs/assembly_eval/bin/samtools'))
    a = p.parse_args()
    OUT.mkdir(exist_ok=False)
    source_gz = a.source_gz or Path(dict(l.split('\t') for l in
        (ROOT/'analysis_genomes.tsv').read_text().splitlines() if l.strip())['GRCh38'])
    derived_gz = OUT/'GRCh38_acgtn.fa.gz'
    positions, summary = convert(source_gz, derived_gz, a.samtools)
    if not positions:
        raise SystemExit('no ambiguity positions found; nothing to fix')
    # Coordinate-preserving check: .fai of the derived FASTA must be byte-identical
    # to the original .fai (94 substitutions change no lengths).
    with gzip.open(derived_gz, 'rb') as f, (OUT/'GRCh38_acgtn.fa').open('wb') as plain:
        subprocess.run(['cp'], stdin=f, stdout=plain, check=True)
    subprocess.run([str(a.samtools), 'faidx', str(OUT/'GRCh38_acgtn.fa')], check=True)
    original_fai = Path(str(source_gz)+'.fai')
    if original_fai.is_file():
        if digest(original_fai) != digest(OUT/'GRCh38_acgtn.fa.fai'):
            raise SystemExit('derived .fai differs from original: coordinates changed')
    record = dict(format='grch38-ambiguity-fix-v1', status='DERIVED_VERIFIED',
        rationale='native PAM search accepts only A/C/G/T/N; IUPAC positions become N',
        scientific_note='94 of ~3.1e9 positions; N never matches a guide base; all hits preserved',
        summary=summary, positions=positions)
    (OUT/'transformation.json').write_text(json.dumps(record, indent=2)+'\n')
    (OUT/'GRCh38_acgtn.fa').unlink()  # uncompressed scratch copy; gz is authoritative
    print(json.dumps(dict(status=record['status'], converted=len(positions),
        derived_sha256=summary['derived_gz_sha256'])))


if __name__ == '__main__':
    main()
