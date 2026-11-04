#!/usr/bin/env python3
"""Check adapter preserves native representatives and rejects unsupported inputs."""
import csv,itertools,os,shutil,subprocess,sys,tempfile
from pathlib import Path
cli,poc,ref,guides,index=sys.argv[1:]
common=[cli,'--search-backend','pam-anchored','--pam-search-bin',poc,'--columba-index',index,'-r',ref,'-G',guides,'-t','1']
for pam in map(''.join,itertools.product('ACGT',repeat=2)):
 for m,b,z,f in [(0,0,2,'.75'),(0,1,2,'.75'),(2,1,1,'.75'),(5,1,2,'.75'),(5,1,0,'.75'),(20,0,0,'0')]:
  native=subprocess.run([poc,index,ref,guides,str(m),str(z),f,pam,str(b)],capture_output=True,text=True,check=True)
  lines=native.stdout.splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'));rows=csv.DictReader(lines[i:],delimiter='\t')
  expected={(r['guide'],r['contig'],r['strand'],r['start'],r['end'],r['cigar']) for r in rows}
  run=subprocess.run(common+['-p',pam,'-m',str(m),'-b',str(b),'-z',str(z),'-f',f],capture_output=True,text=True,check=True)
  got=[]
  for line in run.stdout.splitlines():
   x=line.split('\t');cig=next(z[5:] for z in x[12:] if z.startswith('cg:Z:'));got.append((x[0],x[5],x[4],x[7],x[8],cig))
  assert len(got)==len(set(got)) and set(got)==expected,(pam,m,b,z,f)
for extra in (['-b','2'],['-z','3'],['--no-filter'],['--allow-noncanonical-pam'],['-t','2'],['-w','40'],['--columba-sam','x'],['-f','nan'],['-p','GN']):
 r=subprocess.run(common+extra,text=True,capture_output=True);assert r.returncode and not r.stdout,(extra,r.stderr)
with tempfile.TemporaryDirectory() as td:
 bad=Path(td)/'bad.fa';bad.write_text(Path(ref).read_text()+'\n')
 r=subprocess.run([cli,'--search-backend','pam-anchored','--pam-search-bin',poc,'--columba-index',index,'-r',str(bad),'-G',guides],capture_output=True,text=True)
 assert r.returncode and not r.stdout and 'does not match' in r.stderr
with tempfile.TemporaryDirectory() as td:
 td=Path(td); prefix=td/'index'
 for src in Path(index).parent.glob(Path(index).name+'.*'):
  shutil.copyfile(src,str(prefix)+str(src)[len(index):])
 component=Path(str(prefix)+'.brt');data=component.read_bytes();component.write_bytes(data[:-1]+bytes([data[-1]^1]))
 cmd=[cli,'--search-backend','pam-anchored','--pam-search-bin',poc,'--columba-index',str(prefix),'-r',ref,'-G',guides]
 fake=td/'python3';fake.write_text('#!/bin/sh\nexit 0\n');fake.chmod(0o755)
 env=dict(os.environ,PATH=str(td)+':'+os.environ.get('PATH',''))
 r=subprocess.run(cmd,capture_output=True,text=True,env=env)
 assert r.returncode and not r.stdout and 'index component' in r.stderr
print('96 PAM/option adapter comparisons and 11 rejection checks passed')
