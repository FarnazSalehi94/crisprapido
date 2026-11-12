#!/usr/bin/env python3
"""Scoring-only NRR domain fixtures: every mismatch count and single-gap position."""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--scorer',type=Path,required=True)
    p.add_argument('--python',required=True)
    p.add_argument('--official-python',required=True)
    p.add_argument('--models',required=True)
    p.add_argument('--upstream',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir()
    sys.path.insert(0,str(a.scorer))
    from adapter import encode_tokens
    guide='ACTGAAATCTGTAAGCAGGC'
    change=dict(zip('ACGT','CGTA'))
    rows=[]
    for suffix in ('GG','AG','GA','AA'):
        for first in 'ACGT':
            for mm in range(7):
                target=''.join(change[x] for x in guide[:mm])+guide[mm:]
                pairs=[(guide,target,'none')]
                if mm <= 4:
                    # RNA gap cannot delete a position already counted as mismatch.
                    for pos in range(mm,20):
                        pairs.append((guide,target[:pos]+'-'+target[pos+1:],'RNA'))
                    for pos in range(21):
                        for base in 'ACGT':
                            pairs.append((guide[:pos]+'-'+guide[pos:],target[:pos]+base+target[pos:],'DNA'))
                for g,t,kind in pairs:
                    rows.append(dict(id=len(rows),guide=g+'NGG',target=t+first+suffix,kind=kind,mismatches=mm))
    inputs=a.output/'inputs.jsonl'
    inputs.write_text(''.join(json.dumps(r)+'\n' for r in rows))
    token_hash=hashlib.sha256(b''.join(struct.pack('<24q',*x) for x in encode_tokens(rows))).hexdigest()
    calls=[ [a.python,a.scorer/'worker.py','--model-root',a.models,'--input',inputs,'--output',a.output/'integrated.jsonl','--metrics',a.output/'integrated-metrics.json'],
           [a.official_python,a.scorer/'official_reference.py','--upstream',a.upstream,'--input',inputs,'--output',a.output/'official.jsonl','--metrics',a.output/'official-metrics.json'] ]
    for i,cmd in enumerate(calls):
        with (a.output/(str(i)+'.stdout')).open('x') as out,(a.output/(str(i)+'.stderr')).open('x') as err:
            subprocess.run(list(map(str,cmd)),stdout=out,stderr=err,check=True,timeout=1200)
    own=[json.loads(s) for s in (a.output/'integrated.jsonl').read_text().splitlines()]
    official=[json.loads(s) for s in (a.output/'official.jsonl').read_text().splitlines()]
    if [r['id'] for r in rows]!=[r['id'] for r in own] or [r['id'] for r in own]!=[r['id'] for r in official]:raise ValueError('identity mismatch')
    errors=[abs(r['score']-s['score']) for r,s in zip(own,official)]
    if not all(0<=e<=1e-6 for e in errors):raise ValueError('numeric disagreement')
    if token_hash!=json.loads((a.output/'official-metrics.json').read_text())['categorical_sha256']:raise ValueError('encoding disagreement')
    (a.output/'summary.json').write_text(json.dumps(dict(status='PASS',rows=len(rows),max_absolute_error=max(errors),categorical_sha256=token_hash,scope='16 NRR triplets; ungapped 0..6 mismatches; single RNA/DNA gaps 0..4 mismatches; every feasible gap position and DNA base; computational not biological validation'),indent=2)+'\n')


if __name__=='__main__':main()
