#!/usr/bin/env python3
"""Score actual prototype CIGARs via Rust and an independent matrix oracle."""
import csv,json,re,subprocess,sys
from pathlib import Path
repo=Path(__file__).resolve().parents[3];sys.path.insert(0,str(repo/'scripts/release_validation'));from audit_emitted import independent_cfd, matrices
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'validation'));from index_check import fasta
ref,guide,output,probe=sys.argv[1:];refs=dict(fasta(ref));guides=dict(fasta(guide));mm,pm=matrices(repo);rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
lines=Path(output).read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'));rows=list(csv.DictReader(lines[i:],delimiter='\t'));queries=[];expected=[];classes=set();strands=set();pams=set()
for r in rows:
 g=guides[r['guide']];s=refs[r['contig']][int(r['start']):int(r['end'])];t=s if r['strand']=='+' else rc(s);cigar=r['cigar'];ops=[(int(n),op) for n,op in re.findall(r'(\d+)([=XID])',cigar)];assert ''.join(str(n)+op for n,op in ops)==cigar
 qidx=tidx=0;aligned_g=[];aligned_t=[]
 for n,op in ops:
  for _ in range(n):
   if op=='I':aligned_g.append(g[qidx]);aligned_t.append('-');qidx+=1
   elif op=='D':aligned_g.append('-');aligned_t.append(t[tidx]);tidx+=1
   else:aligned_g.append(g[qidx]);aligned_t.append(t[tidx]);qidx+=1;tidx+=1
 assert qidx==len(g) and tidx==len(t)
 pam=r['pam'];score=independent_cfd(''.join(aligned_g),''.join(aligned_t),ops,pam,mm,pm);expected.append(score);queries.append('\t'.join((g,t,cigar,pam[1:])));strands.add(r['strand']);pams.add(pam[1:]);classes.update(op for _,op in ops)
run=subprocess.run([probe],input='\n'.join(queries)+'\n',capture_output=True,text=True,check=True,cwd=repo);actual=run.stdout.splitlines();assert len(actual)==len(expected)
for i,(a,b) in enumerate(zip(actual,expected)):
 assert (a=='NA' and b is None) or (a!='NA' and b is not None and abs(float(a)-b)<1e-9),(i,a,b,rows[i])
print(json.dumps({'records':len(rows),'scored':sum(x is not None for x in expected),'unscorable':sum(x is None for x in expected),'strands':sorted(strands),'operations':sorted(classes),'pams':sorted(pams)}))
