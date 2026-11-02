#!/usr/bin/env python3
"""Independent bounded, indexed check of mismatch and one/two-base gap settings."""
import argparse, collections, json, random, subprocess
from pathlib import Path
import re
from audit_emitted import parse_row

p=argparse.ArgumentParser(description=__doc__)
for x in ('output','binary','columba','builder','loader','library-path','matrices'): p.add_argument('--'+x,required=True)
a=p.parse_args(); out=Path(a.output).resolve();out.mkdir(exist_ok=False)
rng=random.Random(782761); dna='ACGT'
guides={'g1':'TCTAGCGTTTAACTTCTTCG','g2':'TGCCACCAATACCCTTTGGA'}
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def mutate(s,n):
 s=list(s)
 for i in (0,3,6,9,12,15,18)[:n]:s[i]=dna[(dna.index(s[i])+1)%4]
 return ''.join(s)
def sites(seq, guide, limit, fraction):
 n=len(seq);L=len(guide)
 for strand in '+-':
  ori=seq if strand=='+' else rc(seq)
  for start in range(n):
   for delta in (-2,-1,0,1,2):
    length=L+delta; end=start+length
    if length<=0 or end+3>n or ori[end+1:end+3]!='GG':continue
    target=ori[start:end]; pam=ori[end:end+3]
    if delta==0: variants=[(guide,target)]
    elif delta<0: variants=[(guide[:i]+guide[i-delta:],target) for i in range(L+delta+1)]
    else: variants=[(guide,target[:i]+target[i+delta:]) for i in range(L+1)]
    if any(sum(x!=y for x,y in zip(x,y))<=limit and sum(x==y for x,y in zip(x,y))/L>=fraction for x,y in variants):
     anchor=end if strand=='+' else n-end-3
     yield (strand,anchor,pam)
seqs={};deliberate=[]
for name,g in guides.items():
 for strand in '+-':
  for case in ('exact','mm5','mm6','mm7','ins2','del2','ins2mm3','del2mm3'):
   target=g
   if case.startswith('mm'):target=mutate(g,int(case[2:]))
   if case.startswith('ins2'):target=(mutate(g,3) if case.endswith('mm3') else g);target=target[:8]+target[10:]
   if case.startswith('del2'):target=(mutate(g,3) if case.endswith('mm3') else g);target=target[:8]+'AC'+target[8:]
   flank=''.join(rng.choice(dna) for _ in range(45));tail=''.join(rng.choice(dna) for _ in range(45));ori=flank+target+'AGG'+tail
   seq=ori if strand=='+' else rc(ori); contig=f'{name}_{case}_{"plus" if strand=="+" else "minus"}';seqs[contig]=seq
   deliberate.append((name,contig,strand,45+len(target) if strand=='+' else len(seq)-45-len(target)-3,case))
(out/'ref.fa').write_text(''.join(f'>{k}\n{v}\n' for k,v in seqs.items()))
(out/'guides.fa').write_text(''.join(f'>{k}\n{v}\n' for k,v in guides.items()))
def run(cmd, **kw):return subprocess.run(cmd,check=True,text=True,timeout=120,**kw)
with (out/'build.log').open('w') as log:run([a.loader,'--library-path',a.library_path,a.builder,'-f',str(out/'ref.fa'),'-r',str(out/'index')],stdout=log,stderr=subprocess.STDOUT)
results=[]
for m,f in ((5,.75),(6,.70),(7,.65)):
 label=f'm{m}';tmp=out/(label+'_tmp');tmp.mkdir()
 import os
 env=dict(os.environ,TMPDIR=str(tmp))
 common=['-r',str(out/'ref.fa'),'-G',str(out/'guides.fa'),'-p','GG','-m',str(m),'-b','1','-z','2','-f',str(f),'-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt')]
 with (out/(label+'.paf')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:run([a.binary]+common+['--columba-bin',a.columba,'--columba-index',str(out/'index'),'--keep-columba-sam'],stdout=stdout,stderr=stderr,env=env)
 if m == 5:
  default_args=['-r',str(out/'ref.fa'),'-G',str(out/'guides.fa'),'-t','1','--mismatch-scores',str(Path(a.matrices)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrices)/'pam_scores.txt'),'--columba-bin',a.columba,'--columba-index',str(out/'index')]
  with (out/'defaults.paf').open('w') as stdout,(out/'defaults.stderr').open('w') as stderr:run([a.binary]+default_args,stdout=stdout,stderr=stderr)
  if sorted((out/'defaults.paf').read_text().splitlines()) != sorted((out/'m5.paf').read_text().splitlines()):raise ValueError('default CLI results differ from explicit frozen options')
 truth={(g,c,st,pos) for g,guide in guides.items() for c,seq in seqs.items() for st,pos,pam in sites(seq,guide,m,f)}
 found=[]
 for line in (out/(label+'.paf')).read_text().splitlines():
  x=line.split('\t');st=x[4];found.append((x[0],x[5],st,int(x[8]) if st=='+' else int(x[7])-3))
  _,tags,ops,_,_=parse_row(line,guides,{k:len(v) for k,v in seqs.items()})
  seq=seqs[x[5]];start,end=int(x[7]),int(x[8]);target=seq[start:end];pam=seq[end:end+3] if st=='+' else rc(seq[start-3:start]);target=target if st=='+' else rc(target)
  assert pam[1:]=='GG'
  qi=ti=matches=mismatches=gaps=0
  for length,op in ops:
   if op=='I': qi+=length;gaps+=1;assert length<=2
   elif op=='D': ti+=length;gaps+=1;assert length<=2
   else:
    for _ in range(length):
     equal=guides[x[0]][qi]==target[ti];assert equal==(op=='=');matches+=equal;mismatches+=not equal;qi+=1;ti+=1
  assert qi==20 and ti==len(target) and mismatches<=m and gaps<=1 and matches/20>=f
  assert tags['nm']==('i',str(mismatches))

 foundset=set(found);missing=sorted(truth-foundset);extra=sorted(foundset-truth)
 maps={}
 for sam in tmp.glob('*/output.sam'):
  head=sam.read_text().splitlines(); pg=next(x for x in head if x.startswith('@PG\t'))
  import re
  edit=int(re.search(r' -e (\d+)(?: |$)',pg).group(1));maps[str(edit)]=sum(not x.startswith('@') and not int(x.split('\t')[1])&4 for x in head)
 target=[(g,c,st,pos,case) for g,c,st,pos,case in deliberate if (g,c,st,pos) in truth]
 res={'m':m,'fraction':f,'edit_bound':m+2,'oracle':len(truth),'reported':len(found),'unique':len(foundset),'recovered':len(truth&foundset),'missed':missing,'unexpected':extra,'duplicate_count':len(found)-len(foundset),'independently_audited_records':len(found),'deliberate_eligible':len(target),'deliberate_recovered':sum((g,c,st,pos) in foundset for g,c,st,pos,case in target),'raw_mapped_by_edit':maps}
 results.append(res);print(json.dumps(res,sort_keys=True),flush=True)
(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
if any(x['missed'] or x['unexpected'] or x['duplicate_count'] or x['deliberate_recovered']!=x['deliberate_eligible'] for x in results):raise SystemExit(1)
