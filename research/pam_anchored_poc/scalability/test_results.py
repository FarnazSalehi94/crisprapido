#!/usr/bin/env python3
"""Check the small published measurement artifacts for internal consistency."""
import json
from pathlib import Path
p=Path(__file__).with_name('results');regions=json.loads((p/'regions.json').read_text());runs=json.loads((p/'runs.json').read_text());summary=json.loads((p/'summary.json').read_text());assert len(regions)==len(summary)==7 and len(runs)==42
for region in regions:
 name=region['name'];rs=[x for x in runs if x['region']==name];s=next(x for x in summary if x['region']==name);assert len(rs)==6 and all(x['status']==0 for x in rs)
 for b in ('columba-wfa2','pam-anchored'):
  r=[x for x in rs if x['backend']==b];assert len(r)==3 and len({x['rows'] for x in r})==1 and s['deterministic'][b]
 audit=json.loads((p/'audits'/f'{name}.json').read_text());assert audit['baseline']==audit['experimental']==audit['both'] and not audit['baseline_only'] and not audit['experimental_only'] and not audit['baseline_duplicates'] and not audit['experimental_duplicates']
 if region['bases']==100000:
  o=audit['oracle'];assert o['baseline_recovered']==o['experimental_recovered']==o['expected'] and not o['baseline_missing'] and not o['experimental_missing'] and not o['baseline_unexpected'] and not o['experimental_unexpected']
assert all(json.loads((p/'n_ablation/index_component_equality.json').read_text()).values())
n=json.loads((p/'n_ablation/runs.json').read_text());assert len(n)==12 and all(x['status']==0 for x in n)
for region in ('natural','shadow'):
 x=[r for r in n if r['region']==region and r['backend']=='pam-anchored'];assert len(x)==3 and all(v['states']==1411449 and v['rows']==2 for v in x)
print('42 workflow runs, seven audits and N ablation results consistent')
