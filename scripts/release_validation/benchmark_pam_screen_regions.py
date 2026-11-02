#!/usr/bin/env python3
"""Benchmark unmodified verifier on bounded indexed 2 Mb windows and one guide."""
import argparse,json,os,re,resource,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for n in ('panel','output','binary','columba','matrices','time-bin'):p.add_argument('--'+n,required=True)
a=p.parse_args();panel=Path(a.panel);out=Path(a.output);out.mkdir(exist_ok=False)
# ELMOD2 is the guide that generated the existing 17.65M genome candidates.
guide='ELMOD2_TCTAGCGTTTAACTTCTTCG';(out/'guide.fa').write_text('>'+guide+'\nTCTAGCGTTTAACTTCTTCG\n')
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
def run(cmd,stdout,stderr):
 with stdout.open('w') as o,stderr.open('w') as e:subprocess.run(cmd,check=True,stdout=o,stderr=e,timeout=90,preexec_fn=limit)
def possible(seq,pos,end,rev):
 if rev:lo=max(0,pos-24);hi=min(end-4,len(seq)-21);pair='CC';j=lo;limit=hi+2
 else:lo=max(18,pos+1);hi=min(end+21,len(seq)-3);pair='GG';j=lo+1;limit=hi+3
 if lo>hi:return False
 while True:
  i=seq.find(pair,j,limit)
  if i<0:return False
  pam=i if rev else i-1
  if lo<=pam<=hi and seq[pam+2 if rev else pam] in 'ACGT':return True
  j=i+1
results=[]
for win in json.loads((panel/'windows.json').read_text()):
 d=out/win['name'];d.mkdir();src=panel/win['name'];lines=(src/'ref.fa').read_text().splitlines();seq=''.join(lines[1:]).upper();refname=lines[0][1:].split()[0]
 for edit in (7,0):
  run([a.time_bin,'-v','-o',str(d/f'columba{edit}.time'),a.columba,'-r',str(src/'index'),'-f',str(out/'guide.fa'),'-a','all','-m','edit','-e',str(edit),'-t','1','-o',str(d/f'e{edit}.sam')],d/f'columba{edit}.log',d/f'columba{edit}.err')
 seen=set();counts={'raw':0,'duplicates':0,'direct_invalid':0,'no_overlapping_pam':0,'retained':0};filter_start=time.monotonic()
 with (d/'baseline.sam').open('w') as base,(d/'overlap.sam').open('w') as filt:
  for edit in (7,0):
   for line in (d/f'e{edit}.sam').open():
    if line.startswith('@'):
     if edit==7:base.write(line);filt.write(line)
     continue
    f=line.split('\t');flag=int(f[1])
    if flag&4:continue
    counts['raw']+=1;pos=int(f[3])-1;rev=bool(flag&16);key=(f[0],f[2],pos,rev,f[5])
    if key in seen:counts['duplicates']+=1;continue
    seen.add(key);assert f[2]==refname
    length=sum(int(v) for v,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]) if op in 'MDN=X');end=pos+length
    direct=(pos>=3 and seq[pos-3:pos-1]=='CC' and seq[pos-1] in 'ACGT') if rev else (end+3<=len(seq) and seq[end+1:end+3]=='GG' and seq[end] in 'ACGT')
    allowed=possible(seq,pos,end,rev);counts['direct_invalid']+=not direct;counts['no_overlapping_pam']+=not allowed
    base.write(line)
    if allowed:filt.write(line);counts['retained']+=1
 counts['unique']=len(seen);counts['filter_seconds']=time.monotonic()-filter_start
 modes={}
 for mode in ('baseline','overlap'):
  cmd=[a.time_bin,'-v','-o',str(d/f'{mode}.time'),a.binary,'-r',str(src/'ref.fa'),'-G',str(out/'guide.fa'),'-p','GG','-m','5','-b','1','-z','2','-f','.75','-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-sam',str(d/f'{mode}.sam')]
  t=time.monotonic();run(cmd,d/f'{mode}.paf',d/f'{mode}.stderr');seconds=time.monotonic()-t
  lines=(d/f'{mode}.paf').read_text().splitlines();keys=[]
  for line in lines:
   f=line.split('\t');keys.append((f[0],f[5],f[4],int(f[8]) if f[4]=='+' else int(f[7])-3))
  if len(keys)!=len(set(keys)):raise ValueError('duplicate PAM anchor')
  times=(d/f'{mode}.time').read_text();rss=re.search(r'Maximum resident set size \(kbytes\): (\d+)',times)
  modes[mode]={'seconds':seconds,'max_rss_kib':int(rss.group(1)),'anchors':len(keys),'keys':set(keys),'records':sorted(lines)}
 before,after=modes['baseline'],modes['overlap'];missing=before['keys']-after['keys'];extra=after['keys']-before['keys'];identical=before['records']==after['records']
 for m in modes.values():del m['keys'];del m['records']
 result={'window':win['name'],'region':win['region'],'counts':counts,'modes':modes,'missing_anchors':len(missing),'extra_anchors':len(extra),'records_identical':identical};results.append(result)
 (out/'summary.partial.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(result),flush=True)
 if missing or extra:raise ValueError('overlap screen changes reported anchors')
(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
