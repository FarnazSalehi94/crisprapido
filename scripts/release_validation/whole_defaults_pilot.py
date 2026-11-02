#!/usr/bin/env python3
"""One resource-capped whole-assembly search, audit, exact and baseline comparison."""
import argparse,csv,hashlib,json,os,re,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('root','assembly','guide','reference','index','old-root','binary','columba','matrices','samtools','time-bin','oracle'):p.add_argument('--'+n,required=True)
a=p.parse_args();root=Path(a.root);d=root/'runs'/a.assembly/a.guide;d.mkdir(parents=True,exist_ok=False);(d/'tmp').mkdir()
expected='cb10e5d6be970bb6ef90aac4f2fe5ab035ae2b4d4ce643f28ec32f3895d7e42d'
if hashlib.sha256(Path(a.binary).read_bytes()).hexdigest()!=expected:raise ValueError('development executable hash mismatch')
guides=root/'inputs'/f'{a.guide}.fa';env=dict(os.environ,TMPDIR=str(d/'tmp'))
def run(cmd,stdout,stderr):
 with (d/stdout).open('x') as o,(d/stderr).open('x') as e:subprocess.run(cmd,check=True,stdout=o,stderr=e,env=env)
cmd=[a.time_bin,'-v','-o',str(d/'search.time'),a.binary,'-r',a.reference,'-G',str(guides),'-p','GG','-m','5','-b','1','-z','2','-f','.75','-t','4','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-bin',a.columba,'--columba-index',a.index,'--keep-columba-sam']
start=time.monotonic();run(cmd,'actual.paf','search.stderr');elapsed=time.monotonic()-start
sams={};seen=set();raw={}
for sam in (d/'tmp').glob('*/output.sam'):
 with sam.open() as stream:
  for line in stream:
   if line.startswith('@PG\t'):
    m=re.search(r'(?:^| )-e (\d+)(?: |\t|$)',line)
    if m:
     edit=int(m.group(1))
     if edit in sams:raise ValueError('duplicate edit bound')
     sams[edit]=sam
    break
if set(sams)!={0,7}:raise ValueError('expected Columba e=0 and e=7 SAM')
for edit in (7,0):
 count=0
 with sams[edit].open() as stream:
  for line in stream:
   if line.startswith('@'):continue
   f=line.split('\t');flag=int(f[1])
   if flag&4:continue
   count+=1;seen.add((f[0],f[2],int(f[3]),bool(flag&16),f[5]))
 raw[str(edit)]=count
run(['python3',str(root/'tools/audit_emitted.py'),'--paf',str(d/'actual.paf'),'--fasta',a.reference,'--guides',str(guides),'--matrices',a.matrices,'--samtools',a.samtools,'--output',str(d/'audit'),'--max-mismatches','5','--max-bulges','1','--max-bulge-size','2','--min-match-fraction','.75'],'audit.stdout','audit.stderr')
if not (d/'audit/PASS').exists():raise ValueError('independent auditor did not pass')
run(['python3',str(root/'tools/compare_genome_exact.py'),'--oracle',a.oracle,'--sam',str(sams[0]),'--fasta',a.reference,'--guides',str(guides),'--samtools',a.samtools,'--output',str(d/'exact.json')],'exact.stdout','exact.stderr')
exact=json.loads((d/'exact.json').read_text())
if exact['missed'] or exact['unexpected']:raise ValueError('independent exact discrepancy')
def anchor_set(paf):
 rows=[]
 with Path(paf).open() as f:
  for line in f:
   x=line.split('\t');rows.append((x[0],x[5],x[4],int(x[8]) if x[4]=='+' else int(x[7])-3))
 if len(rows)!=len(set(rows)):raise ValueError('duplicate PAM anchor in '+str(paf))
 return set(rows)
new=anchor_set(d/'actual.paf');old=anchor_set(Path(a.old_root)/'runs'/a.assembly/a.guide/'actual.paf')
truth=[json.loads(x) for x in Path(a.oracle).read_text().splitlines()];truthkeys={(x['guide'],x['contig'],x['strand'],x['end'] if x['strand']=='+' else x['start']-3) for x in truth}
if truthkeys-new:raise ValueError('independent exact anchor absent from final output')
comparison={'old':len(old),'new':len(new),'shared':len(old&new),'new_only':len(new-old),'old_only':len(old-new),'old_only_keys':sorted(old-new)}
(d/'new_only.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in sorted(new-old)))
summary={'assembly':a.assembly,'guide':a.guide,'binary_sha256':expected,'search_wall_seconds':elapsed,'raw_mapped_e7':raw['7'],'raw_mapped_e0':raw['0'],'unique_candidate_keys':len(seen),'reported_anchors':len(new),'wfa2_calls':None,'wfa2_calls_reason':'uninstrumented development executable; replay would repeat costly verification','comparison':comparison,'exact':exact}
(d/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
if old-new:raise ValueError('previous eligible PAM anchors absent; investigate before continuing')
(d/'PASS').write_text('PASS\n');print(json.dumps(summary,sort_keys=True))
