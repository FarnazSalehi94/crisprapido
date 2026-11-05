#!/usr/bin/env python3
"""Audit integrated backend on the four independently enumerated N substitutions."""
import argparse,json,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--panel',required=True);p.add_argument('--crispr',required=True);p.add_argument('--poc',required=True);p.add_argument('--out',required=True);a=p.parse_args();panel=Path(a.panel);out=Path(a.out);out.mkdir(parents=True,exist_ok=True);root=Path(__file__).resolve().parents[3];auditor=root/'research/pam_anchored_poc/integration/audit.py';results={};prior=json.loads((panel/'summary.json').read_text())
for c in 'ACGT':
 idx=panel/c/'index';paf=out/(c+'.paf');cmd=[a.crispr,'--search-backend','pam-anchored','--pam-search-bin',a.poc,'--columba-index',str(idx),'-r',str(panel/'ref.fa'),'-G',str(panel/'guides.fa'),'-t','1','-m','5','-b','1','-z','2','-f','.75','-p','GG']
 with paf.open('w') as f,(out/(c+'.log')).open('w') as e:subprocess.run(cmd,stdout=f,stderr=e,check=True,timeout=60,cwd=root)
 audit=subprocess.run([sys.executable,str(auditor),str(panel/'ref.fa'),str(panel/'guides.fa'),str(panel/'oracle.jsonl'),str(paf)],capture_output=True,text=True,check=True,cwd=root);results[c]={'integrated':json.loads(audit.stdout),'baseline':prior['replacements'][c]['baseline'],'baseline_misses_with_candidate_on_same_contig':prior['replacements'][c]['baseline_misses_with_candidate_on_same_contig']}
(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,sort_keys=True))
