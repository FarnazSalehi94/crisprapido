#!/usr/bin/env python3
"""Differential test of all four FM-index substitutes for reference N."""
import argparse,collections,json,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'validation'))
from index_check import check
DNA='ACGT';rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def change(s,positions):
 s=list(s)
 for p in positions:s[p]=DNA[(DNA.index(s[p])+1)%4]
 return ''.join(s)
def main():
 a=argparse.ArgumentParser()
 for x in ('out','builder','poc','baseline','columba','matrix-dir'):a.add_argument('--'+x,required=True)
 a=a.parse_args();out=Path(a.out);out.mkdir(exist_ok=False,parents=True)
 guides={'complex':'ACGTTGCACTGAGTACCGTA','repeat':'AC'*10};refs={};design=[]
 for gid,g in guides.items():
  # Five mismatches count the aligned N itself; insertion cases use <=4
  # because f=.75 requires 15 matched guide bases.
  base=change(g,[1,4,12,16]); base=base[:8]+'N'+base[9:]
  cases={'n_mm5':base,'n_mm5_D1':base[:10]+'T'+base[10:], 'n_mm5_D2':base[:10]+'TA'+base[10:],
         'n_mm4_I1':(change(g,[1,12,16])[:8]+'N'+change(g,[1,12,16])[9:])[:5]+(change(g,[1,12,16])[:8]+'N'+change(g,[1,12,16])[9:])[6:],
         'n_mm3_I2':(change(g,[1,16])[:8]+'N'+change(g,[1,16])[9:])[:5]+(change(g,[1,16])[:8]+'N'+change(g,[1,16])[9:])[7:],
         'three_N':g[:7]+'NNN'+g[10:], 'two_N_in_D2':g[:9]+'NN'+g[9:]}
  for case,target in cases.items():
   for strand in '+-':
    name=gid+'_'+case+('_plus' if strand=='+' else '_minus');oriented='CC'+target+'AGGG'+'ACACAC';genomic=oriented if strand=='+' else rc(oriented);refs[name]=genomic
    anchor=(2+len(target)) if strand=='+' else len(genomic)-(2+len(target))-3
    design.append({'guide':gid,'contig':name,'strand':strand,'pam_start':anchor,'class':case})
 (out/'ref.fa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in refs.items()));(out/'guides.fa').write_text(''.join('>'+k+'\n'+s+'\n' for k,s in guides.items()))
 truth={}
 for c,seq in refs.items():
  for strand in '+-':
   ori=seq if strand=='+' else rc(seq)
   for end in range(18,len(ori)-2):
    if ori[end] not in DNA or ori[end+1:end+3]!='GG':continue
    anchor=end if strand=='+' else len(seq)-end-3
    for delta in (-2,-1,0,1,2):
     length=20+delta
     if end<length:continue
     t=ori[end-length:end]
     for gid,g in guides.items():
      if delta==0:variations=[(g,t,'M')]
      elif delta<0:variations=((g[:i]+g[i-delta:],t,'I'+str(-delta)) for i in range(length+1))
      else:variations=((g,t[:i]+t[i+delta:],'D'+str(delta)) for i in range(21))
      for q,r,op in variations:
       mismatches=sum(x!=y for x,y in zip(q,r));matches=len(q)-mismatches
       if mismatches<=5 and matches>=15:
        truth.setdefault((gid,c,strand,anchor),set()).add((op,mismatches,'N' in t))
 missing=[x for x in design if (x['guide'],x['contig'],x['strand'],x['pam_start']) not in truth]
 assert not missing,missing
 required={'n_mm5':('M',5,True),'n_mm5_D1':('D1',5,True),'n_mm5_D2':('D2',5,True),'n_mm4_I1':('I1',4,True),'n_mm3_I2':('I2',3,True),'three_N':('M',3,True),'two_N_in_D2':('D2',0,True)}
 absent=[x for x in design if required[x['class']] not in truth[(x['guide'],x['contig'],x['strand'],x['pam_start'])]];assert not absent,absent
 (out/'oracle.jsonl').write_text(''.join(json.dumps({'guide':g,'contig':c,'strand':s,'anchor':pos,'categories':sorted(op+'_'+str(mm)+'_N'+str(int(n)) for op,mm,n in v)},sort_keys=True)+'\n' for (g,c,s,pos),v in sorted(truth.items())))
 (out/'design.json').write_text(json.dumps(design,indent=2)+'\n')
 summary={}
 for letter in DNA:
  sub=out/letter;sub.mkdir();prefix=sub/'index';idxfa=sub/'index_ref.fa';idxfa.write_text(''.join('>'+k+'\n'+s.replace('N',letter)+'\n' for k,s in refs.items()))
  subprocess.run([a.builder,'-r',str(prefix),'-f',str(idxfa)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
  (sub/'index-check.json').write_text(json.dumps(check(prefix,out/'ref.fa'),indent=2)+'\n')
  for name in ('ref.fa','guides.fa','oracle.jsonl'):(sub/name).symlink_to((out/name).resolve())
  with (sub/'poc.tsv').open('w') as f,(sub/'poc.log').open('w') as e:subprocess.run([a.poc,str(prefix),str(out/'ref.fa'),str(out/'guides.fa'),'5','2','.75'],stdout=f,stderr=e,check=True)
  cmd=[a.baseline,'-r',str(out/'ref.fa'),'-G',str(out/'guides.fa'),'-p','GG','-m','5','-b','1','-z','2','-f','.75','-t','1','--mismatch-scores',str(Path(a.matrix_dir)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrix_dir)/'pam_scores.txt'),'--columba-bin',a.columba,'--columba-index',str(prefix)]
  with (sub/'baseline.paf').open('w') as f,(sub/'baseline.log').open('w') as e:subprocess.run(cmd,stdout=f,stderr=e,check=True,timeout=120)
  result=subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'validation'/'compare.py'),str(sub)],capture_output=True,text=True,check=True)
  data=json.loads(result.stdout);(sub/'comparison.json').write_text(json.dumps(data,indent=2)+'\n')
  assert data['prototype']['missed']==0 and data['prototype']['unexpected']==0 and data['prototype']['duplicates']==0, data['prototype']
  summary[letter]={k:{v:data[k][v] for v in ('expected','reported','recovered','missed','unexpected','duplicates')} for k in ('prototype','baseline')}
  # Distinguish discovery from downstream verification for missed contigs.
  found=set()
  for distance in (7,0):
   sam=sub/('raw_e'+str(distance)+'.sam')
   subprocess.run([a.columba,'-r',str(prefix),'-f',str(out/'guides.fa'),'-a','all','-m','edit','-e',str(distance),'-t','1','-o',str(sam)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
   for line in sam.read_text().splitlines():
    if line.startswith('@'):continue
    f=line.split('\t')
    if not int(f[1])&4: found.add((f[0],f[2]))
  missing=data['baseline']['missing_keys']
  summary[letter]['baseline_misses_with_candidate_on_same_contig']=sum((g,c) in found for g,c,s,anchor in missing)
  detail=[]
  for gid,contig,strand,anchor in missing:
   seq=refs[contig].replace('N',letter);ori=seq if strand=='+' else rc(seq);end=anchor if strand=='+' else len(seq)-anchor-3;g=guides[gid];best=100
   for delta in (-2,-1,0,1,2):
    length=20+delta
    if end<length:continue
    target=ori[end-length:end]
    if delta==0:variants=[(g,target)]
    elif delta<0:variants=((g[:i]+g[i-delta:],target) for i in range(length+1))
    else:variants=((g,target[:i]+target[i+delta:]) for i in range(21))
    for q,t in variants:best=min(best,abs(delta)+sum(x!=y for x,y in zip(q,t)))
   detail.append({'key':[gid,contig,strand,anchor],'indexed_one_gap_edit_distance':best,'has_any_candidate_on_contig':(gid,contig) in found})
   assert best<=7
  summary[letter]['baseline_missing_details']=detail
 print(json.dumps({'cases':len(design),'classes':list(required),'oracle_anchors':len(truth),'replacements':summary},sort_keys=True))
 (out/'summary.json').write_text(json.dumps({'cases':len(design),'oracle_anchors':len(truth),'replacements':summary},indent=2)+'\n')
if __name__=='__main__':main()
