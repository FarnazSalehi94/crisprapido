#!/usr/bin/env python3
"""Make deterministic 1 kb windows from the first two indexed contigs per assembly."""
import argparse
import csv
import json
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--manifest', required=True)
p.add_argument('--samtools', required=True)
p.add_argument('--output', required=True)
a = p.parse_args()
out = Path(a.output)
out.mkdir(exist_ok=False)
records = []
sequences = []
with open(a.manifest) as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        fasta = row['fasta']
        fai = Path(fasta + '.fai')
        for i, line in enumerate(fai.open()):
            if i == 2:
                break
            contig, length = line.split('\t')[:2]
            if int(length) < 101000:
                raise ValueError('contig shorter than fixed interval')
            region = contig + ':100001-101000'
            data = subprocess.run([a.samtools, 'faidx', fasta, region], check=True,
                                  stdout=subprocess.PIPE, universal_newlines=True).stdout
            sequence = ''.join(data.splitlines()[1:]).upper()
            if len(sequence) != 1000:
                raise ValueError('incomplete interval: ' + region)
            name = row['assembly_id'] + '_window_' + str(i + 1)
            sequences.append('>' + name + '\n' + sequence + '\n')
            records.append(dict(window=name, source_fasta=fasta, region=region,
                                source_start_zero_based=100000, length=len(sequence)))
if len(records) != 24:
    raise ValueError('expected twelve assemblies with two contigs each')
(out / 'ref.fa').write_text(''.join(sequences))
(out / 'windows.json').write_text(json.dumps(records, indent=2) + '\n')
