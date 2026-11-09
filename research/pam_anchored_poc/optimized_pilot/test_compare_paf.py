#!/usr/bin/env python3
"""Exercise exact set comparison on the retained, already audited pilot output."""
import json,subprocess,sys,tempfile
from pathlib import Path
baseline=Path('/export/local/home/salehi/crisprapido-resource-pilot-results/hg00735-pilot-2871278/output.paf');tool=Path(__file__).with_name('compare_paf.py');lines=baseline.read_text().splitlines();assert len(lines)==21447
with tempfile.TemporaryDirectory() as tmp:
 tmp=Path(tmp);cand=tmp/'candidate.paf';out=tmp/'result.json'
 def run(rows,success):
  cand.write_text('\n'.join(rows)+'\n');p=subprocess.run([sys.executable,str(tool),'--baseline',str(baseline),'--candidate',str(cand),'--out',str(out)],capture_output=True,text=True);assert (p.returncode==0)==success,(p.stdout,p.stderr);return json.loads(out.read_text())
 x=run(lines,True);assert x['shared']==21447 and x['representative_alignment_differences']==0 and x['cfd_differences']==0
 x=run(lines+[lines[0]],False);assert x['candidate_duplicate_anchors']==1
 x=run(lines[:-1],False);assert x['missing']==1
 altered=lines.copy();i=next(i for i,x in enumerate(altered) if '\tcf:f:' in x);fields=altered[i].split('\t');j=next(j for j,x in enumerate(fields) if x.startswith('cf:f:'));fields[j]='cf:f:9.9999';altered[i]='\t'.join(fields);x=run(altered,True);assert x['cfd_differences']==1
print('21447-site comparison detects duplicates, missing anchors and CFD differences')
