#!/usr/bin/env python3
"""Check compact results for internal consistency without genomic inputs."""
import json
from pathlib import Path
p=Path(__file__).with_name('results');read=lambda x:json.loads((p/x).read_text())
r=read('3mb/runs.json');assert len(r)==4 and all(x['status']==0 and x['stop_reason'] is None for x in r)
for region,count in [('CHM13',20),('GRCh38',18)]:
 a=read('3mb/'+region+'-audit.json');assert a['baseline']==a['experimental']==a['both']==count and not a['baseline_only'] and not a['experimental_only'] and not a['baseline_duplicates'] and not a['experimental_duplicates']
 assert len([x for x in r if x['region']==region and x['rows']==count])==2
n=read('n_substitution/integrated.json')
for b in 'ACGT':
 x=n[b]['integrated'];assert x['expected']==x['reported']==x['recovered']==40 and x['missing']==x['unexpected']==x['duplicates']==0
v=read('verify-comparison.json');assert sum(x['parser']=='old' for x in v)==sum(x['parser']=='new' for x in v)==5
assert read('172-guide-characteristics.json')['N_overlapping_sites']==75
print('3 Mb, N-substitution, verifier and 172-panel results consistent')
