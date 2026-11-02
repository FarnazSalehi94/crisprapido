#!/usr/bin/env python3
"""One bounded Slurm pilot: frozen indexed search, diagnostic replay, independent audit."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
for name in ('root', 'assembly', 'guide', 'reference', 'index', 'binary', 'profile', 'columba', 'matrices', 'samtools', 'time-bin', 'oracle'):
    p.add_argument('--' + name, required=True)
a = p.parse_args()
root = Path(a.root).resolve()
if not (root / 'verification/PASS').is_file():
    raise RuntimeError('full manifest verification has not passed')
if hashlib.sha256(Path(a.binary).read_bytes()).hexdigest() != 'c37f4ded6e7c6afb3a7d22eccc5566f291244731bf6872dbfbd0b4d47f53aa74':
    raise RuntimeError('frozen binary hash mismatch')
if hashlib.sha256(Path(a.profile).read_bytes()).hexdigest() != 'fef24de98613802c344b8dddff6f776e60827677bc548df27c733619c813057a':
    raise RuntimeError('diagnostic binary hash mismatch')
freeze = json.loads((root/'inputs/freeze_manifest.json').read_text())
for path, key in [(a.columba, 'columba_wrapper_sha256'), (str(Path(a.columba).with_name('columba')), 'columba_sha256'), (a.samtools, 'samtools_sha256'), (str(Path(a.matrices)/'mismatch_scores.txt'), 'mismatch_matrix_sha256'), (str(Path(a.matrices)/'pam_scores.txt'), 'pam_matrix_sha256')]:
    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != freeze[key]:
        raise RuntimeError('frozen dependency hash mismatch: ' + path)
d = root / 'runs' / a.assembly / a.guide
d.mkdir(parents=True, exist_ok=False)
(d / 'tmp').mkdir()
guides = root / 'inputs' / (a.guide + '.fa')
env = dict(os.environ, TMPDIR=str(d / 'tmp'))
common = ['-r', a.reference, '-G', str(guides), '-p', 'GG', '-m', '4', '-b', '1', '-z', '1', '-f', '0.75', '-t', '4', '--mismatch-scores', str(Path(a.matrices)/'mismatch_scores.txt'), '--pam-scores', str(Path(a.matrices)/'pam_scores.txt')]
def run(cmd, out, err):
    with (d/out).open('x') as stdout, (d/err).open('x') as stderr:
        subprocess.run(cmd, stdout=stdout, stderr=stderr, env=env, check=True)
start = time.monotonic()
run([a.time_bin, '-v', '-o', str(d/'search.time'), a.binary] + common + ['--columba-bin', a.columba, '--columba-index', a.index, '--keep-columba-sam'], 'actual.paf', 'search.stderr')
search_wall = time.monotonic() - start
# Preserve first-seen ordering and exactly match the application's key and unmapped handling.
sams = {}
for sam in (d/'tmp').glob('*/output.sam'):
    with sam.open() as f:
        for line in f:
            if line.startswith('@PG\t'):
                match = re.search(r'(?:^| )-e (\d+)(?: |\t|$)', line)
                if match:
                    edit = int(match.group(1))
                    if edit in sams: raise ValueError('duplicate SAM edit bound')
                    sams[edit] = sam
                break
if set(sams) != {0, 5}: raise ValueError('expected one e=5 and one e=0 SAM: '+repr(sams))
seen = set(); counts = {}
with (d/'union.sam').open('x') as out:
    for edit in (5, 0):
        count = 0
        with sams[edit].open() as f:
            for line in f:
                if line.startswith('@'):
                    if edit == 5: out.write(line)
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 11: raise ValueError('malformed SAM row')
                flag = int(fields[1])
                if flag & 4: continue
                count += 1
                key = (fields[0], fields[2], int(fields[3]), bool(flag & 16), fields[5])
                if key not in seen:
                    seen.add(key); out.write(line)
        counts[str(edit)] = count
run([a.profile] + common + ['--columba-sam', str(d/'union.sam')], 'profile.paf', 'profile.stderr')
# The program can emit records in a different order with four threads.
def records(path):
    with path.open() as f: return sorted(line.rstrip('\n') for line in f if line.strip())
actual, replay = records(d/'actual.paf'), records(d/'profile.paf')
if actual != replay: raise ValueError('instrumented replay differs from frozen output')
profile = {}
with (d/'profile.stderr').open() as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) == 4 and fields[0] == 'PROFILE': profile[fields[1]] = {'calls':int(fields[2]), 'nanoseconds':int(fields[3])}
if 'wfa_align' not in profile: raise ValueError('missing WFA2 counter')
tools = root/'tools'
run(['python3', str(tools/'audit_emitted.py'), '--paf', str(d/'actual.paf'), '--fasta', a.reference, '--guides', str(guides), '--matrices', a.matrices, '--samtools', a.samtools, '--output', str(d/'audit')], 'audit.stdout', 'audit.stderr')
if not (d/'audit/PASS').exists(): raise ValueError('audit did not pass')
run(['python3', str(tools/'compare_genome_exact.py'), '--oracle', a.oracle, '--sam', str(sams[0]), '--fasta', a.reference, '--guides', str(guides), '--samtools', a.samtools, '--output', str(d/'exact.json')], 'exact.stdout', 'exact.stderr')
exact = json.loads((d/'exact.json').read_text())
if exact['missed'] or exact['unexpected']: raise ValueError('exact comparison discrepancy')
# Require exact-truth PAM anchors in the final audited output as well as raw e=0.
truth = [json.loads(line) for line in Path(a.oracle).read_text().splitlines()]
truth_anchors = {(x['guide'], x['contig'], x['strand'], x['end'] if x['strand']=='+' else x['start']-3) for x in truth}
with (d/'audit/records.tsv').open() as f:
    audited = list(csv.DictReader(f, delimiter='\t'))
anchors = {(x['guide'], x['contig'], x['strand'], int(x['end0']) if x['strand']=='+' else int(x['start0'])-3) for x in audited}
missing_anchors = truth_anchors - anchors
if missing_anchors: raise ValueError('exact truth anchors absent from final output: '+repr(sorted(missing_anchors)))
if len(anchors) != len(audited): raise ValueError('duplicate emitted PAM anchors')
exact['final_exact_truth_anchors_recovered'] = len(truth_anchors & anchors)
summary = {'assembly':a.assembly, 'guide':a.guide, 'search_wall_seconds':search_wall, 'raw_mapped_e5':counts['5'], 'raw_mapped_e0':counts['0'], 'unique_candidate_keys':len(seen), 'wfa2_calls_replay':profile['wfa_align']['calls'], 'reported_records':len(actual), 'exact':exact, 'profile':profile, 'search_time_file':str(d/'search.time')}
(d/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
(d/'PASS').write_text('PASS\n')
print(json.dumps(summary, sort_keys=True))
