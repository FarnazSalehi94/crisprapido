#!/usr/bin/env python3
"""Validate adjacent I/D event accounting against production validation semantics."""
import csv,re,subprocess,sys,json
from pathlib import Path
builder,poc,out=sys.argv[1:];p=Path(out);p.mkdir(exist_ok=False);g='ACGTTGCACTGAGTACCGTA';t=g[:9]+'A'+g[10:];assert t[9]!=g[9]
(p/'ref.fa').write_text('>plus\n'+t+'AGG\n>minus\n'+'CCT'+t.translate(str.maketrans('ACGT','TGCA'))[::-1]+'\n');(p/'guide.fa').write_text('>g\n'+g+'\n')
# Match the current fail-closed native index-binding contract.
subprocess.run([sys.executable,str(Path(__file__).with_name('index_manifest.py')),'build','--builder',builder,'--reference',str(p/'ref.fa'),'--index',str(p/'index')],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
r=subprocess.run([poc,str(p/'index'),str(p/'ref.fa'),str(p/'guide.fa'),'0','2','.75','GG','1'],text=True,capture_output=True,check=True);(p/'poc.tsv').write_text(r.stdout);l=r.stdout.splitlines();i=next(i for i,x in enumerate(l) if x.startswith('guide\tcontig\t'));rows=list(csv.DictReader(l[i:],delimiter='\t'))
for c,s,a in [('plus','+',20),('minus','-',3)]:
 row=next(x for x in rows if x['contig']==c and x['strand']==s and int(x['anchor'])==a)
 ops=re.findall(r'(\d+)([=XID])',row['cigar']);assert int(row['mismatches'])==0 and int(row['bulge_events'])==1 and int(row['inserted_guide_bases'])==1 and int(row['deleted_guide_bases'])==1 and int(row['match_count'])==19
 assert any(ops[i][1] in 'ID' and ops[i+1][1] in 'ID' and ops[i][1]!=ops[i+1][1] for i in range(len(ops)-1))
print(json.dumps({'adjacent_mixed_gap_strands':2,'mismatches':0,'events':1,'combined_gap_length':2,'all_pass':True}))
