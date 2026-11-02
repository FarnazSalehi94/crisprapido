#!/usr/bin/env python3
"""Paired timing and exact-output comparison of disposable interval-PAM probe."""
import argparse,json,re,resource,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('panel','regions','real-panel','benchmark','output','baseline','probe','matrices','time-bin'):p.add_argument('--'+n,required=True)
a=p.parse_args();out=Path(a.output);out.mkdir(exist_ok=False)
def limit():resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
cases=[('oracle_8kb',Path(a.real_panel)/'ref.fa',Path(a.real_panel)/'guides.fa',Path(a.benchmark)/'baseline.sam',Path(a.benchmark)/'baseline.paf')]
for w in json.loads((Path(a.panel)/'windows.json').read_text()):
 n=w['name'];cases.append((n,Path(a.panel)/n/'ref.fa',Path(a.regions)/'guide.fa',Path(a.regions)/n/'baseline.sam',Path(a.regions)/n/'baseline.paf'))
results=[]
for name,ref,guides,sam,expected_path in cases:
 d=out/name;d.mkdir();expected=sorted(expected_path.read_text().splitlines());modes={}
 for mode,exe in (('baseline',a.baseline),('interval_probe',a.probe)):
  cmd=[a.time_bin,'-v','-o',str(d/f'{mode}.time'),exe,'-r',str(ref),'-G',str(guides),'-p','GG','-m','5','-b','1','-z','2','-f','.75','-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-sam',str(sam)]
  start=time.monotonic()
  with (d/f'{mode}.paf').open('w') as o,(d/f'{mode}.stderr').open('w') as e:subprocess.run(cmd,check=True,stdout=o,stderr=e,preexec_fn=limit,timeout=90)
  elapsed=time.monotonic()-start;observed=sorted((d/f'{mode}.paf').read_text().splitlines());identical=observed==expected
  txt=(d/f'{mode}.time').read_text();rss=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)',txt).group(1));modes[mode]={'seconds':elapsed,'max_rss_kib':rss,'rows':len(observed),'identical_to_reference':identical}
  if not identical:raise ValueError(name+' '+mode+' output differs')
 row={'case':name,'modes':modes};results.append(row);(out/'summary.partial.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(row),flush=True)
(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
