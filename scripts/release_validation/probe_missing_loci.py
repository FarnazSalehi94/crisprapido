#!/usr/bin/env python3
"""Replay each missing PAM-anchored oracle interval as a separate imported candidate."""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path
from resolve_bulge_discrepancies import analyze, fasta, locus
from coverage_oracle import enumerate_sites

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('directory')
p.add_argument('--binary', required=True)
p.add_argument('--output', required=True)
a = p.parse_args()
d = Path(a.directory).resolve()
refs, guides = fasta(d/'ref.fa'), fasta(d/'guides.fa')
_, metrics = analyze(d)
missing = {tuple(x) for x in metrics['missing_loci']}
keys = set()
for guide, sequence in guides.items():
    for contig, ref in refs.items():
        for row in enumerate_sites(ref, sequence):
            key = (guide, contig, row['start'], row['end'], row['strand'])
            if locus(key, refs) in missing:
                keys.add(key)
records = []
with tempfile.TemporaryDirectory() as temp:
    sam = Path(temp)/'candidate.sam'
    for key in sorted(keys):
        guide, contig, start, end, strand = key
        # Only reference span and strand are used by this anchored diagnostic;
        # this deliberately supplied record is not a Columba discovery result.
        sam.write_text('guide\t{}\t{}\t{}\t60\t{}M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0\n'.format(
            16 if strand == '-' else 0, contig, start+1, end-start))
        base = [a.binary, '-r', str(d/'ref.fa'), '-g', guides[guide], '--columba-sam', str(sam),
                '-p','GG','-m','4','-b','1','-z','1','-f','0.75','-t','1']
        attempts = {}
        for label, extra in [('strict', []), ('allow_noncanonical', ['--allow-noncanonical-pam']),
                             ('no_filter', ['--no-filter'])]:
            result = subprocess.run(base+extra, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, check=True)
            attempts[label] = {'paf':result.stdout.splitlines(), 'stderr':result.stderr.splitlines()}
        records.append({'key':key, 'pam_locus':locus(key, refs), 'attempts':attempts})
Path(a.output).write_text(json.dumps(records, indent=2)+'\n')
print(json.dumps({'intervals':len(records), 'strict_with_paf':sum(bool(x['attempts']['strict']['paf']) for x in records)}))
