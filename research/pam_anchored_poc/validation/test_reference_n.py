#!/usr/bin/env python3
"""Independent bounded exhaustive regression for reference N and PAM boundaries."""
import csv,json,random,re,subprocess,sys
from pathlib import Path
R=random.Random(17293);dna='ACGT';rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
builder,poc,out=sys.argv[1:];p=Path(out);p.mkdir(exist_ok=False)
g='ACGTTGCACTGAGTACCGTA';assert len(g)==20
with_n=lambda t,i:t[:i]+'N'+t[i+1:]
cases={
 'target_n':(with_n(g,10),'AGG',True),
 'n_near_pam':(with_n(g,19),'TGG',True),
 'multiple_n':(g[:9]+'NNN'+g[12:],'CGG',True),
 'n_near_insertion':(with_n(g[:7]+g[9:],8),'GGG',True),
 'n_near_deletion':(with_n(g[:8]+'AC'+g[8:],10),'AGG',True),
 'n_inside_deletion':(g[:8]+'NN'+g[8:],'TGG',True),
 'pam_wildcard_N':(g,'NGG',False),
 'pam_first_G_N':(g,'ANG',False),
 'pam_second_G_N':(g,'AGN',False),
}
refs={}
for name,(target,pam,_) in cases.items():
 refs[name+'_plus']=target+pam
 refs[name+'_minus']=rc(target+pam)
(p/'ref.fa').write_text(''.join('>'+n+'\n'+s+'\n' for n,s in refs.items()));(p/'guides.fa').write_text('>g\n'+g+'\n')
subprocess.run([builder,'-r',str(p/'index'),'-f',str(p/'ref.fa')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
with (p/'poc.tsv').open('w') as o,(p/'poc.log').open('w') as e:subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(p/'guides.fa'),'5','2','.75'],check=True,stdout=o,stderr=e)
truth={}
for c,seq in refs.items():
 for strand in '+-':
  ori=seq if strand=='+' else rc(seq)
  for end in range(18,len(ori)-2):
   pam=ori[end:end+3]
   if pam[0] not in dna or pam[1:]!='GG':continue
   for delta in (-2,-1,0,1,2):
    n=20+delta
    if end<n:continue
    t=ori[end-n:end]
    if delta==0: variants=[(g,t,'ungapped')]
    elif delta<0:variants=((g[:i]+g[i-delta:],t,'I'+str(-delta)) for i in range(n+1))
    else:variants=((g,t[:i]+t[i+delta:],'D'+str(delta)) for i in range(21))
    for q,v,cat in variants:
     mm=sum(a!=b for a,b in zip(q,v))
     if mm<=5 and len(q)-mm>=15:
      key=('g',c,strand,end if strand=='+' else len(seq)-end-3);truth.setdefault(key,set()).add(cat)
ls=(p/'poc.tsv').read_text().splitlines();i=next(i for i,x in enumerate(ls) if x.startswith('guide\tcontig\t'));rows=list(csv.DictReader(ls[i:],delimiter='\t'));got={}
for x in rows:
 start,end=int(x['start']),int(x['end']);seq=refs[x['contig']];target=seq[start:end] if x['strand']=='+' else rc(seq[start:end]);qi=ti=mm=match=ins=dele=0
 ops=re.findall(r'(\d+)([=XID])',x['cigar']);assert ''.join(n+op for n,op in ops)==x['cigar']; gaps=0
 for n,op in ops:
  n=int(n)
  if op in 'ID':gaps+=1;assert n<=2
  for _ in range(n):
   if op=='I':qi+=1;ins+=1
   elif op=='D':ti+=1;dele+=1
   else:
    assert qi<len(g) and ti<len(target);same=g[qi]==target[ti];assert same==(op=='=');qi+=1;ti+=1;mm+=not same;match+=same
 assert qi==20 and ti==len(target) and mm<=5 and match>=15 and gaps<=1
 assert (mm,ins,dele)==tuple(int(x[k]) for k in ('mismatches','inserted_guide_bases','deleted_guide_bases'))
 key=(x['guide'],x['contig'],x['strand'],int(x['anchor'])-(3 if x['strand']=='-' else 0)); assert key not in got;got[key]=x
assert set(got)==set(truth),(list(set(truth)-set(got)),list(set(got)-set(truth)))
for name,(target,pam,valid) in cases.items():
 for suffix,strand in (('plus','+'),('minus','-')):
  key=('g',name+'_'+suffix,strand,len(target) if strand=='+' else 0)
  assert (key in truth)==valid,(name,strand,key)
  if valid:
   x=got[key]; assert x['pam']==pam and x['cigar'] and int(x['mismatches'])<=5
   if name in ('target_n','n_near_pam'):assert int(x['mismatches'])==1
   if name=='multiple_n':assert int(x['mismatches'])==3
   if name=='n_inside_deletion':assert int(x['mismatches'])==0 and int(x['deleted_guide_bases'])==2
   if name=='n_near_insertion':assert 'I2' in truth[key]
   if name in ('n_near_deletion','n_inside_deletion'):assert 'D2' in truth[key]
# At m=0 a reference N can be consumed by an adjacent I/D event or by a
# reference-only D; neither consumes a mismatch. The no-gap test below
# independently checks that an aligned N consumes the mismatch budget.
with (p/'zero.tsv').open('w') as o,(p/'zero.log').open('w') as e:
 subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(p/'guides.fa'),'0','2','.75'],check=True,stdout=o,stderr=e)
lines=(p/'zero.tsv').read_text().splitlines(); j=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'));zero=list(csv.DictReader(lines[j:],delimiter='\t'))
for strand in '+-':
 suffix='plus' if strand=='+' else 'minus'
 assert any(x['contig']=='target_n_'+suffix and x['strand']==strand and int(x['inserted_guide_bases'])==1 and int(x['deleted_guide_bases'])==1 for x in zero)
 assert any(x['contig']=='n_inside_deletion_'+suffix and x['strand']==strand and int(x['mismatches'])==0 and int(x['deleted_guide_bases'])==2 for x in zero)
# Force a replacement-equals-guide case: use the indexed replacement at the
# first target's N as the guide base. True replay must still charge mismatch.
indexed=(p/'index.txt.bin').read_bytes()[4:];offset=cases['target_n'][0].index('N');replacement=chr(indexed[offset]);assert replacement in dna
masked_guide=g[:offset]+replacement+g[offset+1:];(p/'replacement_guide.fa').write_text('>replacement\n'+masked_guide+'\n')
for budget in (0,1):
 with (p/('replacement'+str(budget)+'.tsv')).open('w') as o,(p/('replacement'+str(budget)+'.log')).open('w') as e:
  subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(p/'replacement_guide.fa'),str(budget),'0','.75'],check=True,stdout=o,stderr=e)
 lines=(p/('replacement'+str(budget)+'.tsv')).read_text().splitlines(); j=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'));found=list(csv.DictReader(lines[j:],delimiter='\t'))
 selected=[x for x in found if x['contig']=='target_n_plus' and x['strand']=='+' and int(x['anchor'])==20]
 if budget==0:assert not selected
 else:assert len(selected)==1 and int(selected[0]['mismatches'])==1 and 'X' in selected[0]['cigar']
print(json.dumps({'cases':len(cases)*2,'expected_anchors':len(truth),'recovered':len(got),'missed':len(set(truth)-set(got)),'unexpected':len(set(got)-set(truth)),'duplicates':len(rows)-len(got),'invalid_pam_cases':6,'strands':['+','-'],'seed':17293},sort_keys=True))
