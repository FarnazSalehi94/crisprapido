#!/usr/bin/env python3
"""Fresh-process wall/CPU/ru_maxrss wrapper. Use the tree supervisor separately."""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

p=argparse.ArgumentParser()
for key in ('metrics','stdout','stderr'):
    p.add_argument('--'+key,type=Path,required=True)
p.add_argument('--stdin',type=Path)
p.add_argument('--timeout',type=float,default=300)
p.add_argument('command',nargs=argparse.REMAINDER)
a=p.parse_args()
if a.command and a.command[0]=='--':a.command=a.command[1:]
if not a.command or not math.isfinite(a.timeout) or a.timeout<=0:raise SystemExit('invalid measurement arguments')
start=time.perf_counter();reason='completed'
with a.stdout.open('xb') as stdout,a.stderr.open('xb') as stderr:
    source=a.stdin.open('rb') if a.stdin else open(os.devnull,'rb')
    try:
        child=subprocess.Popen(a.command,stdin=source,stdout=stdout,stderr=stderr,start_new_session=True)
        try:code=child.wait(timeout=a.timeout)
        except subprocess.TimeoutExpired:
            reason='timeout';os.killpg(child.pid,signal.SIGKILL);code=child.wait()
    finally:source.close()
usage=resource.getrusage(resource.RUSAGE_CHILDREN)
with a.metrics.open('x') as f:
    json.dump(dict(reason=reason,returncode=code,wall_seconds=time.perf_counter()-start,
        child_user_seconds=usage.ru_utime,child_system_seconds=usage.ru_stime,
        max_individual_child_rss_kib=usage.ru_maxrss,
        memory_scope='Linux max individual descendant RSS, not simultaneous process-tree memory'),f,indent=2)
    f.write('\n')
if code or reason!='completed':sys.exit(1)
