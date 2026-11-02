#!/usr/bin/env python3
"""Audit disposable interval probe output on bounded reference windows."""
import argparse,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('panel','regions','paired','output','samtools','matrices','scripts'):p.add_argument('--'+n,required=True)
a=p.parse_args();out=Path(a.output);out.mkdir(exist_ok=False);results=[]
for w in json.loads((Path(a.panel)/'windows.json').read_text()):
 n=w['name'];d=out/n;d.mkdir();ref=d/'ref.fa';ref.symlink_to((Path(a.panel)/n/'ref.fa').resolve());subprocess.run([a.samtools,'faidx',str(ref)],check=True,timeout=60)
 cmd=['python3',str(Path(a.scripts)/'audit_emitted.py'),'--paf',str(Path(a.paired)/n/'interval_probe.paf'),'--fasta',str(ref),'--guides',str(Path(a.regions)/'guide.fa'),'--matrices',a.matrices,'--samtools',a.samtools,'--output',str(d/'audit'),'--max-mismatches','5','--max-bulges','1','--max-bulge-size','2','--min-match-fraction','.75']
 with (d/'stdout').open('w') as o,(d/'stderr').open('w') as e:subprocess.run(cmd,check=True,stdout=o,stderr=e,timeout=120)
 summary=json.loads((d/'audit/summary.json').read_text());results.append({'window':n,'status':summary['status'],'rows':summary['rows'],'errors':summary['errors']})
(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
