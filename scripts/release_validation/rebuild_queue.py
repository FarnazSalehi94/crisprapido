#!/usr/bin/env python3
"""Bounded queue for explicitly listed failed index sets, never genome queries."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','evidence-dir','builder','release-manifest','output']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--assemblies',nargs='+',required=True)
    p.add_argument('--existing-status',nargs='*',default=[])
    p.add_argument('--max-builders',type=int,default=2)
    p.add_argument('--start-after-status')
    a=p.parse_args()
    if not 1<=a.max_builders<=2: p.error('max-builders must be 1 or 2')
    out=Path(a.output);out.mkdir(exist_ok=False)
    if a.start_after_status:
        while not Path(a.start_after_status).exists():time.sleep(5)
    waiting=list(a.assemblies);active=[];finished=[]
    while waiting or active:
        existing=0
        for path in a.existing_status:
            try:
                existing+=json.loads(Path(path).read_text())['status'] in ['PREPARING','BUILDING']
            except (OSError, ValueError):
                existing+=1  # incomplete status is conservatively occupied
        for item in list(active):
            if item['process'].poll() is not None:
                item['log'].close();finished.append(dict(assembly=item['assembly'],exit_code=item['process'].returncode,command=item['command']))
                active.remove(item)
                (out/'finished.json').write_text(json.dumps(finished,indent=2)+'\n')
        while waiting and len(active)+existing<a.max_builders:
            asm=waiting.pop(0)
            command=[sys.executable,str(Path(__file__).with_name('rebuild_index.py').resolve()),'--evidence',str(Path(a.evidence_dir)/(asm+'.json')),'--assembly',asm,'--builder',a.builder,'--release-manifest',a.release_manifest,'--output',str(Path(a.root)/('rebuilt_'+asm))]
            log=(out/(asm+'.log')).open('x')
            process=subprocess.Popen(command,stdout=log,stderr=log)
            active.append(dict(assembly=asm,process=process,log=log,command=command))
            (out/(asm+'.command.json')).write_text(json.dumps(command,indent=2)+'\n')
            print('started '+asm,flush=True)
        if waiting or active:time.sleep(5)
    ok=all(item['exit_code']==0 for item in finished)
    (out/'summary.json').write_text(json.dumps(dict(status='PASS' if ok else 'FAIL',finished=finished,max_builders=a.max_builders),indent=2)+'\n')
    return 0 if ok else 1


if __name__=='__main__': sys.exit(main())
