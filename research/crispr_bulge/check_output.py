#!/usr/bin/env python3
"""Independently compare annotated PAF scores with saved official predictions.

No model execution, reference loading, CFD recalculation or search. This checks
full postprocessor output in addition to direct worker prediction parity.
"""
import argparse
import json
import math
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--validation-dir',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
root=a.validation_dir
summary=json.loads((root/'summary.json').read_text())
expected=[json.loads(line) for line in (root/'official.jsonl').read_text().splitlines()]
assert [x['id'] for x in expected]==list(range(len(expected)))
result={}
for mode in ('both','crispr-bulge'):
    scores=[];rows=0
    for panel in summary['panels']:
        directory=root/(panel['name']+'-'+mode)
        meta=json.loads((directory/'metadata.json').read_text())
        text=(directory/'output.paf').read_text().splitlines()
        assert len(text)==panel['rows']
        assert meta['statuses']==panel['statuses']
        for line in text:
            f=line.split('\t');tags={x.split(':',2)[0]:x.split(':',2)[2] for x in f[12:]}
            assert ('cb' in tags)==(tags['bu']=='scored')
            if mode=='crispr-bulge':assert 'cf' not in tags
            if 'cb' in tags:
                value=float(tags['cb']);assert math.isfinite(value) and 0<=value<=1
                scores.append(value)
        rows+=len(text)
    assert len(scores)==len(expected)
    maximum=max(abs(score-reference['score']) for score,reference in zip(scores,expected))
    assert maximum<=1e-6,maximum
    result[mode]=dict(input_rows_preserved=rows,predictions_compared=len(scores),
        unavailable_rows_preserved=rows-len(scores),maximum_absolute_difference=maximum,tolerance=1e-6)
with a.output.open('x') as f:
    json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps(result,sort_keys=True))
