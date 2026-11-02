#!/usr/bin/env python3
"""Bounded real-human indexed comparison against direct PAM-anchored enumeration."""
import argparse, collections, hashlib, json, os, random, re, subprocess, time
from pathlib import Path
from audit_emitted import parse_row

p=argparse.ArgumentParser(description=__doc__)
for n in ('source-panel','output','binary','columba','builder','loader','library-path','samtools','matrices','time-bin'):p.add_argument('--'+n,required=True)
a=p.parse_args();source=Path(a.source_panel);out=Path(a.output).resolve();out.mkdir(exist_ok=False)
DNA='ACGT'; trans=str.maketrans('ACGTN','TGCAN');rc=lambda s:s.translate(trans)[::-1]
def run(cmd,**kw):return subprocess.run(cmd,check=True,text=True,timeout=300,**kw)
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
 return h.hexdigest()
windows=json.loads((source/'windows.json').read_text());seqs={};meta=[]
for w in windows:
 text=(source/w['name']/'ref.fa').read_text();lines=text.splitlines();seq=''.join(lines[1:]).upper();assert len(seq)==2000000
 # Measure within-parent-window repeat proxy with nonoverlapping 20mers, then
 # choose one 2kb block; this is not a formal repeat annotation.
 counts=collections.Counter(seq[i:i+20] for i in range(0,len(seq)-19,20))
 blocks=[]
 for off in range(0,len(seq)-2000+1,2000):
  block=seq[off:off+2000];repeat=sum(counts[block[i:i+20]]>1 for i in range(0,1981,20))
  if set(block)<=set(DNA):blocks.append((repeat,off))
 repeat,off=(max(blocks) if 'repeat_candidate' in w['name'] and 'nonrepeat' not in w['name'] else min(blocks))
 chrom,span=w['region'].rsplit(':',1);begin=int(span.split('-')[0])+off;region=f'{chrom}:{begin}-{begin+1999}'
 assembly=w['assembly']; refs={'CHM13':'/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/recovered_text_CHM13/input/chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz','GRCh38':'/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/rebuilt_GRCh38/input/GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz'}
 fetched=run([a.samtools,'faidx',refs[assembly],region],stdout=subprocess.PIPE).stdout
 fresh=''.join(fetched.splitlines()[1:]).upper();assert fresh==seq[off:off+2000]
 seqs[w['name']]=fresh;meta.append({'name':w['name'],'assembly':assembly,'reference':refs[assembly],'region':region,'offset_0based':begin-1,'sampled_repeated_20mers':repeat,'sampled_20mers':100,'sha256':hashlib.sha256(fresh.encode()).hexdigest()})
(out/'ref.fa').write_text(''.join(f'>{k}\n{v}\n' for k,v in seqs.items()))
# Derive 20-nt guides independently from reference targets at genuine PAMs.
guides={};design=[]
def change(s,n):
 z=list(s)
 for i in [1,4,7,11,15][:n]:z[i]=DNA[(DNA.index(z[i])+1)%4]
 return ''.join(z)
for contig,seq in seqs.items():
 for strand in '+-':
  orient=seq if strand=='+' else rc(seq)
  anchors=[e for e in range(50,len(seq)-50) if orient[e+1:e+3]=='GG']
  if not anchors:raise ValueError('no PAM for '+contig)
  e=anchors[len(anchors)//2];pam=orient[e:e+3]
  for case in ('exact','mm1','mm2','mm3','mm4','mm5','I1','I2','D1','D2','I2mm2','D2mm2'):
   if case.startswith('mm'):g=change(orient[e-20:e],int(case[2:]))
   elif case=='exact':g=orient[e-20:e]
   elif case.startswith('I'):
    gap=int(case[1]);t=orient[e-(20-gap):e];g=t[:8]+'AC'[:gap]+t[8:]
    if 'mm2' in case:g=change(g,2)
   else:
    gap=int(case[1]);t=orient[e-(20+gap):e];g=t[:8]+t[8+gap:]
    if 'mm2' in case:g=change(g,2)
   assert len(g)==20 and set(g)<=set(DNA)
   name=f'{contig}_{"plus" if strand=="+" else "minus"}_{case}';guides[name]=g
   anchor=e if strand=='+' else len(seq)-e-3
   design.append({'guide':name,'contig':contig,'strand':strand,'anchor':anchor,'case':case,'pam':pam})
(out/'guides.fa').write_text(''.join(f'>{k}\n{v}\n' for k,v in guides.items()))
(out/'panel.json').write_text(json.dumps({'windows':meta,'designs':design},indent=2)+'\n')
with (out/'build.log').open('w') as log:run([a.loader,'--library-path',a.library_path,a.builder,'-f',str(out/'ref.fa'),'-r',str(out/'index')],stdout=log,stderr=subprocess.STDOUT)
# Independent exhaustive enumeration of every PAM, span and possible one-gap
# placement, comparing letters directly (no alignment library or Columba).
t0=time.monotonic();truth={};categories=collections.Counter()
for contig,seq in seqs.items():
 for strand in '+-':
  ori=seq if strand=='+' else rc(seq)
  pams=[e for e in range(18,len(ori)-2) if ori[e+1:e+3]=='GG']
  for e in pams:
   for delta in (-2,-1,0,1,2):
    length=20+delta
    if e<length:continue
    target=ori[e-length:e]
    if set(target)-set(DNA):continue
    anchor=e if strand=='+' else len(seq)-e-3
    for name,g in guides.items():
     if delta==0:variants=[(g,target)]
     elif delta<0:variants=((g[:i]+g[i-delta:],target) for i in range(length+1))
     else:variants=((g,target[:i]+target[i+delta:]) for i in range(21))
     for q,t in variants:
      mm=sum(a!=b for a,b in zip(q,t));matches=len(q)-mm
      if mm<=5 and matches>=15:
       key=(name,contig,strand,anchor);cat=('exact' if mm==0 and delta==0 else 'mismatch_'+str(mm) if delta==0 else ('I' if delta<0 else 'D')+str(abs(delta))+('_mm'+str(mm) if mm else ''))
       truth.setdefault(key,set()).add(cat)
# Every design is expected to be an eligible site, with intended structure.
missing_design=[d for d in design if (d['guide'],d['contig'],d['strand'],d['anchor']) not in truth]
if missing_design:raise ValueError('ineligible designed site: '+repr(missing_design))
def intended_category(case):
 if case.startswith('mm'):return 'mismatch_'+case[2:]
 if 'mm' in case:return case.split('mm')[0]+'_mm'+case.split('mm')[1]
 return case
missing_intended=[d for d in design if intended_category(d['case']) not in truth[(d['guide'],d['contig'],d['strand'],d['anchor'])]]
if missing_intended:raise ValueError('intended structure not independently eligible: '+repr(missing_intended))
for cats in truth.values():categories.update(cats)
(out/'oracle.jsonl').write_text(''.join(json.dumps({'guide':g,'contig':c,'strand':s,'anchor':pos,'categories':sorted(v)},sort_keys=True)+'\n' for (g,c,s,pos),v in sorted(truth.items())))
oracle_seconds=time.monotonic()-t0
cmd=[a.time_bin,'-v','-o',str(out/'search.time'),a.binary,'-r',str(out/'ref.fa'),'-G',str(out/'guides.fa'),'-p','GG','-m','5','-b','1','-z','2','-f','0.75','-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-bin',a.columba,'--columba-index',str(out/'index'),'--keep-columba-sam']
tmp=out/'tmp';tmp.mkdir();env=dict(os.environ,TMPDIR=str(tmp));t0=time.monotonic()
with (out/'actual.paf').open('w') as stdout,(out/'search.stderr').open('w') as stderr:run(cmd,stdout=stdout,stderr=stderr,env=env)
search_seconds=time.monotonic()-t0
found=[];audit=[]
for line in (out/'actual.paf').read_text().splitlines():
 x=line.split('\t');strand=x[4];start,end=int(x[7]),int(x[8]);key=(x[0],x[5],strand,end if strand=='+' else start-3);found.append(key)
 _,tags,ops,_,_=parse_row(line,guides,{k:len(v) for k,v in seqs.items()})
 seq=seqs[x[5]];target=seq[start:end];pam=seq[end:end+3] if strand=='+' else rc(seq[start-3:start]);target=target if strand=='+' else rc(target)
 if len(pam)!=3 or pam[1:]!='GG':raise ValueError('invalid PAM '+line)
 qi=ti=matches=mm=gaps=0
 for length,op in ops:
  if op=='I':qi+=length;gaps+=1;assert length<=2
  elif op=='D':ti+=length;gaps+=1;assert length<=2
  else:
   for _ in range(length):
    same=guides[x[0]][qi]==target[ti];assert same==(op=='=');matches+=same;mm+=not same;qi+=1;ti+=1
 assert qi==20 and ti==len(target) and gaps<=1 and mm<=5 and matches>=15 and tags['nm']==('i',str(mm))
 audit.append({'key':key,'cigar':tags['cg'][1],'mismatches':mm,'gaps':gaps,'pam':pam})
foundset=set(found);truthset=set(truth);missing=sorted(truthset-foundset);unexpected=sorted(foundset-truthset)
counts={}
for cat in sorted(categories):
 keys={k for k,v in truth.items() if cat in v};counts[cat]={'expected':len(keys),'recovered':len(keys&foundset)}
groups={}
for contig in seqs:
 for strand in '+-':
  keys={k for k in truthset if k[1]==contig and k[2]==strand}
  groups[contig+'_'+strand]={'expected':len(keys),'recovered':len(keys&foundset)}
res={'implementation_commit':'49f20fc','binary_sha256':digest(a.binary),'windows':len(seqs),'total_bp':sum(map(len,seqs.values())),'guides':len(guides),'expected':len(truth),'reported':len(found),'unique_reported':len(foundset),'recovered':len(truthset&foundset),'recall':len(truthset&foundset)/len(truthset),'precision':len(truthset&foundset)/len(foundset) if foundset else None,'duplicates':len(found)-len(foundset),'missed':missing,'unexpected':unexpected,'designed_total':len(design),'designed_recovered':sum((d['guide'],d['contig'],d['strand'],d['anchor']) in foundset for d in design),'categories':counts,'groups':groups,'oracle_seconds':oracle_seconds,'search_seconds':search_seconds,'independently_audited':len(audit)}
(out/'audit.jsonl').write_text(''.join(json.dumps(v,sort_keys=True)+'\n' for v in audit));(out/'summary.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,sort_keys=True))
if missing or unexpected or len(found)!=len(foundset):raise SystemExit(1)
