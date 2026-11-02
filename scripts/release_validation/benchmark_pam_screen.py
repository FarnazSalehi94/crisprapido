#!/usr/bin/env python3
"""Compare experimental PAM candidate screens with a bounded independent oracle.

All modes use the same unmodified development binary. Imported candidate SAMs
are experimental inputs; discovery and filtering times are reported separately.
"""
import argparse,bisect,collections,json,os,re,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('panel','output','binary','columba','matrices','time-bin','samtools','scripts'):p.add_argument('--'+n,required=True)
a=p.parse_args();panel=Path(a.panel);out=Path(a.output);out.mkdir(exist_ok=False)
def run(cmd,**kw):return subprocess.run(cmd,check=True,text=True,timeout=300,**kw)
refs={};name=''
for line in (panel/'ref.fa').read_text().splitlines():
 if line.startswith('>'):name=line[1:];refs[name]=''
 else:refs[name]+=line.strip().upper()
pams={}
for name,seq in refs.items():
 pams[name,'+']=[i for i in range(len(seq)-2) if seq[i] in 'ACGT' and seq[i+1:i+3]=='GG']
 pams[name,'-']=[i for i in range(len(seq)-2) if seq[i:i+2]=='CC' and seq[i+2] in 'ACGT']
def contains(name,strand,lo,hi):
 arr=pams[name,strand];i=bisect.bisect_left(arr,lo);return i<len(arr) and arr[i]<=hi
samlist=list((panel/'tmp').glob('*/output.sam'));sams={}
for path in samlist:
 for line in path.open():
  if line.startswith('@PG\t'):
   m=re.search(r'(?:^| )-e (\d+)(?: |\t|$)',line)
   if m:sams[int(m.group(1))]=path
   break
if set(sams)!={0,7}:raise ValueError('need real indexed e=0 and e=7 SAMs')
# Compare a deliberately tightened raw edit bound as a counterexample, not a
# proposed substitute: it has no general completeness guarantee.
t0=time.monotonic()
with (out/'e6.log').open('w') as o,(out/'e6.stderr').open('w') as e:
 run([a.time_bin,'-v','-o',str(out/'e6.time'),a.columba,'-r',str(panel/'index'),'-f',str(panel/'guides.fa'),'-a','all','-m','edit','-e','6','-t','1','-o',str(out/'e6.raw.sam')],stdout=o,stderr=e)
discovery_e6=time.monotonic()-t0
truth={(x['guide'],x['contig'],x['strand'],x['anchor']) for x in map(json.loads,(panel/'oracle.jsonl').read_text().splitlines())}
def prepare(paths,mode):
 seen=set();raw=0;removed=0;by=collections.Counter();path=out/(mode+'.sam');t=time.monotonic()
 if any(path.resolve()==src.resolve() for src in paths):raise ValueError('refusing to overwrite an input SAM')
 with path.open('x') as dst:
  for index,src in enumerate(paths):
   for line in src.open():
    if line.startswith('@'):
     if index==0:dst.write(line)
     continue
    f=line.split('\t');flag=int(f[1]);
    if flag&4:continue
    raw+=1;pos=int(f[3])-1;rev=bool(flag&16);key=(f[0],f[2],pos,rev,f[5])
    if key in seen:removed+=1;continue
    seen.add(key);span=sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]) if op in 'MDN=X');end=pos+span;seq=refs[f[2]];strand='-' if rev else '+'
    direct=contains(f[2],strand,pos-3 if rev else end,pos-3 if rev else end)
    if rev:possible=contains(f[2],strand,max(0,pos-24),min(end-4,len(seq)-21))
    else:possible=contains(f[2],strand,max(18,pos+1),min(end+21,len(seq)-3))
    by['direct_pam_invalid']+=not direct;by['no_overlapping_pam']+=not possible
    if mode=='direct_pam' and not direct:continue
    if mode=='overlap_pam' and not possible:continue
    dst.write(line);by['retained']+=1
 return {'raw_mapped':raw,'unique':len(seen),'duplicates':removed,'screen_counts':dict(by),'seconds':time.monotonic()-t}
results={}
for mode,paths in [('baseline',[sams[7],sams[0]]),('direct_pam',[sams[7],sams[0]]),('overlap_pam',[sams[7],sams[0]]),('e6',[out/'e6.raw.sam',sams[0]])]:
 stats=prepare(paths,mode);cmd=[a.time_bin,'-v','-o',str(out/(mode+'.time')),a.binary,'-r',str(panel/'ref.fa'),'-G',str(panel/'guides.fa'),'-p','GG','-m','5','-b','1','-z','2','-f','.75','-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-sam',str(out/(mode+'.sam'))]
 t=time.monotonic()
 with (out/(mode+'.paf')).open('w') as stdout,(out/(mode+'.stderr')).open('w') as stderr:run(cmd,stdout=stdout,stderr=stderr)
 stats['verify_seconds']=time.monotonic()-t;found=[]
 for line in (out/(mode+'.paf')).read_text().splitlines():
  f=line.split('\t');found.append((f[0],f[5],f[4],int(f[8]) if f[4]=='+' else int(f[7])-3))
 keys=set(found);miss=truth-keys;extra=keys-truth
 stats.update(oracle_expected=len(truth),reported=len(found),unique_reported=len(keys),recovered=len(keys&truth),recall=len(keys&truth)/len(truth),precision=len(keys&truth)/len(keys) if keys else None,missed=[list(k) for k in sorted(miss)],unexpected=[list(k) for k in sorted(extra)],duplicate_output=len(found)-len(keys))
 with (out/(mode+'.audit.stdout')).open('w') as o,(out/(mode+'.audit.stderr')).open('w') as e:run(['python3',str(Path(a.scripts)/'audit_emitted.py'),'--paf',str(out/(mode+'.paf')),'--fasta',str(panel/'ref.fa'),'--guides',str(panel/'guides.fa'),'--matrices',a.matrices,'--samtools',a.samtools,'--output',str(out/(mode+'.audit')),'--max-mismatches','5','--max-bulges','1','--max-bulge-size','2','--min-match-fraction','.75'],stdout=o,stderr=e)
 stats['audit']=json.loads((out/(mode+'.audit/summary.json')).read_text())['status'];results[mode]=stats
 print(json.dumps({'mode':mode,**stats},sort_keys=True),flush=True)
results['e6_discovery_seconds']=discovery_e6
(out/'summary.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
