#!/usr/bin/env python3
"""Scoring-only saved pilot validation. No search command exists here."""
import collections
import csv
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys


def main():
    root = Path(sys.argv[1]).resolve()
    p = json.loads((root/'plan.json').read_text())
    sys.path.insert(0, str(Path(p['worker']).parent))
    import worker as w
    w.verify_pins(p)
    sys.path.insert(0, str(root/'scorer'))
    from adapter import encode_tokens
    results = {}
    requests = []
    token_hash = hashlib.sha256()
    for pam in ('NGG','NAG','NGA','NAA'):
        marker = json.loads((Path(p['pilot_root'])/'01_pilot'/pam/'SUCCESS.json').read_text())
        for f,h in marker['outputs'].items():
            if w.digest(f) != h: raise ValueError('pilot output changed: '+f)
        raw = next(f for f in marker['outputs'] if f.endswith('/search.stdout'))
        old = next(f for f in marker['outputs'] if f.endswith('/scoring.stdout'))
        audited = next(f for f in marker['outputs'] if f.endswith('/audit/records.tsv'))
        out = root/pam; out.mkdir()
        command = [p['cli'], '-r', p['reference'], '-G', p['guides'], '--score-paf', raw,
          '--scoring','both','--scoring-output-dir',out/'neural', '--crispr-bulge-python',p['model_python'],
          '--crispr-bulge-model-root',p['models'],'--scoring-script',root/'scorer/score_paf.py',
          '--scoring-reference-sha256',p['reference_sha256'],'--scoring-reference-index',p['reference_index'],
          '--scoring-reference-index-sha256',p['pins'][p['reference_index']],'--scoring-worker-timeout',1800]
        w.run(command,out/'both.paf',out/'stderr',out/'metrics.json',p,root)
        # Exercise the unchanged CLI default scoring path on every saved input.
        w.run([p['cli'],'-r',p['reference'],'-G',p['guides'],'--score-paf',raw],out/'default.paf',out/'default.stderr',out/'default.metrics.json',p,root)
        if w.digest(raw) != w.digest(out/'default.paf'): raise ValueError('default CFD drift')
        if pam == 'NGG' and w.digest(old) != w.digest(out/'both.paf'): raise ValueError('NGG scoring drift')
        counts = collections.Counter(); strata = collections.Counter()
        with open(raw) as original, (out/'both.paf').open() as annotated:
            for a,b in itertools.zip_longest(original,annotated):
                if a is None or b is None or not b.startswith(a.rstrip('\n')+'\t'): raise ValueError('original PAF not preserved')
        with open(audited) as af, (out/'neural/scoring/annotations.tsv').open() as nf:
            for a,n in itertools.zip_longest(csv.DictReader(af,delimiter='\t'),csv.DictReader(nf,delimiter='\t')):
                if a is None or n is None: raise ValueError('row loss')
                if a['line'] != str(int(n['row_index'])+1): raise ValueError('row identity')
                for k in ['contig','strand','pam','aligned_guide','aligned_target']:
                    if a[k] != n[k]: raise ValueError('independent input mismatch: '+k)
                if a['guide'] != n['guide_id'] or a['start0'] != n['start_0based'] or a['end0'] != n['end_0based']: raise ValueError('coordinate mismatch')
                counts['rows'] += 1
                counts['cfd_defined'] += n['cfd'] != 'N/A'
                counts[n['crispr_bulge_reason']] += 1
                if n['crispr_bulge_reason'] != 'scored':
                    if n['crispr_bulge'] != 'N/A': raise ValueError('unavailable numeric value')
                    continue
                row=dict(id=len(requests),guide=a['aligned_guide']+'NGG',target=a['aligned_target']+a['pam'],expected=float(n['crispr_bulge']),pam=pam)
                token_hash.update(struct.pack('<24q',*encode_tokens([row])[0]))
                requests.append(row)
                gap='DNA' if '-' in a['aligned_guide'] else 'RNA' if '-' in a['aligned_target'] else 'none'
                strata[str((a['mismatches'],gap,a['strand']))]+=1
        results[pam]=dict(counts=dict(counts),strata=dict(strata),metrics=json.loads((out/'metrics.json').read_text()))
    with (root/'requests.jsonl').open('x') as f:
        for r in requests:f.write(json.dumps(r)+'\n')
    w.run([p['official_python'],root/'scorer/official_reference.py','--upstream',p['upstream'],'--input',root/'requests.jsonl','--output',root/'official.jsonl','--metrics',root/'official-metrics.json'],root/'official.stdout',root/'official.stderr',root/'official-run-metrics.json',p,root)
    predictions=[json.loads(s) for s in (root/'official.jsonl').read_text().splitlines()]
    if [r['id'] for r in requests] != [r['id'] for r in predictions]:raise ValueError('official identities')
    errors=[abs(r['expected']-s['score']) for r,s in zip(requests,predictions)]
    if not errors or any(not math.isfinite(x) or x>1e-6 for x in errors):raise ValueError('official score disagreement')
    metrics=json.loads((root/'official-metrics.json').read_text())
    if token_hash.hexdigest()!=metrics['categorical_sha256']:raise ValueError('official features disagree')
    w.verify_pins(p)
    w.save(root/'summary.json',dict(status='PASS',runs=results,official_rows=len(requests),max_absolute_error=max(errors),categorical_sha256=token_hash.hexdigest(),limitation='Computational agreement only; noncanonical biological validity not established',outputs={str(f):w.digest(f) for f in root.rglob('*') if f.is_file() and f.name not in ['monitor.log','monitor.json'] and not f.name.startswith('slurm-')}))


if __name__ == '__main__':main()
