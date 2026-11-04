#!/usr/bin/env python3
"""Repeat bounded same-input runs; do not interpret phase timings as identical work."""
import argparse,csv,json,os,re,signal,subprocess,time
from pathlib import Path
from index_check import check
p=argparse.ArgumentParser();
for name in ('reference','index','guides','poc','baseline','columba','matrix_dir','time_bin','out'):p.add_argument('--'+name.replace('_','-'),required=True)
p.add_argument('--keep-columba-logs',action='store_true');p.add_argument('--repetitions',type=int,default=2);p.add_argument('--timeout',type=int,default=120);a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
(out/'index-check.json').write_text(json.dumps(check(a.index,a.reference),indent=2)+'\n')
def keys(path,kind):
 lines=path.read_text().splitlines()
 if kind=='poc':
  i=next(i for i,x in enumerate(lines) if x.startswith('guide\tcontig\t'));return [(x['guide'],x['contig'],x['strand'],int(x['anchor'])) for x in csv.DictReader(lines[i:],delimiter='\t')]
 out=[]
 for l in lines:
  f=l.split('\t');out.append((f[0],f[5],f[4],int(f[8] if f[4]=='+' else f[7])))
 return out
res=[];sets={}
for i in range(a.repetitions):
 for kind in ('poc','baseline'):
  outfile=out/(kind+str(i)+('.tsv' if kind=='poc' else '.paf'));err=out/(kind+str(i)+'.log');tf=out/(kind+str(i)+'.time')
  if kind=='poc': cmd=[a.poc,a.index,a.reference,a.guides,'5','2','.75']
  else: cmd=[a.baseline,'-r',a.reference,'-G',a.guides,'-p','GG','-m','5','-b','1','-z','2','-f','.75','-t','1','--mismatch-scores',str(Path(a.matrix_dir)/'mismatch_scores.txt'),'--pam-scores',str(Path(a.matrix_dir)/'pam_scores.txt'),'--columba-bin',a.columba,'--columba-index',a.index]
  env=os.environ.copy(); tmp=out/('columba_tmp'+str(i));
  if kind=='baseline' and a.keep_columba_logs:
   tmp.mkdir();env['TMPDIR']=str(tmp.resolve());cmd.append('--keep-columba-sam')
  with outfile.open('w') as fo,err.open('w') as fe:
   child=subprocess.Popen([a.time_bin,'-f','%e %U %S %M','-o',str(tf)]+cmd,stdout=fo,stderr=fe,start_new_session=True,env=env)
   try:status=child.wait(timeout=a.timeout)
   except subprocess.TimeoutExpired:
    os.killpg(child.pid,signal.SIGKILL);child.wait();status='timeout'
  if status!=0:res.append({'kind':kind,'run':i,'status':status});continue
  wall,user,system,rss=tf.read_text().split(); lines=err.read_text().splitlines();item={'kind':kind,'run':i,'status':status,'wall_seconds':float(wall),'user_seconds':float(user),'system_seconds':float(system),'peak_process_rss_kib':int(rss)}
  if kind=='poc':
   m=re.search(r'states=(\d+) search_wall_seconds=([\deE.+-]+)',err.read_text()); assert m;item.update(states=int(m[1]),search_wall_seconds=float(m[2])); replay=re.search(r'locate_replay_ns=(\d+)',err.read_text());
   if replay:item.update(locate_replay_seconds=int(replay[1])/1e9, search_excluding_locate_replay_seconds=max(0.0,float(m[2])-int(replay[1])/1e9))
  else:
   item['profile']={x[1]:{'calls':int(x[2]),'wall_ns':int(x[3])} for line in lines if line.startswith('PROFILE\t') for x in [line.split('\t')]}
   if a.keep_columba_logs:
    logs=list(tmp.glob('*/columba.log'));sam=list(tmp.glob('*/output.sam'));item['sam_bytes']=sum(x.stat().st_size for x in sam)
    worker=[float(v) for f in logs for v in re.findall(r'mapping threads finished mapping in ([0-9.]+)s',f.read_text())];item['columba_mapping_worker_seconds']=sum(worker) if len(worker)==len(logs) and logs else None
    item['columba_processes']=len(logs)
  k=keys(outfile,kind);item.update(rows=len(k),unique=len(set(k)),duplicates=len(k)-len(set(k)));sets[(kind,i)]=set(k);res.append(item)
comparison=[]
for i in range(a.repetitions):
 if ('poc',i) in sets and ('baseline',i) in sets:
  x,y=sets['poc',i],sets['baseline',i];comparison.append({'run':i,'poc_only':[list(z) for z in sorted(x-y)],'baseline_only':[list(z) for z in sorted(y-x)]})
(out/'summary.json').write_text(json.dumps({'runs':res,'anchor_comparison':comparison},indent=2)+'\n');print(json.dumps({'runs':res,'anchor_comparison':comparison},sort_keys=True))
