#!/usr/bin/env python3
"""Validate one search per strand per guide and summarize its time and states."""
import json,re,sys
from pathlib import Path
log=Path(sys.argv[1]).read_text();audit=json.loads(Path(sys.argv[2]).read_text());metrics=json.loads(Path(sys.argv[3]).read_text());expected=set(audit['per_guide']);found={};pattern=r'^GUIDE_PROFILE\t([^\t]+)\t([+-])\t(\d+)\t([0-9.eE+-]+)$'
for name,strand,states,seconds in re.findall(pattern,log,re.M):
 assert name in expected and (name,strand) not in found;found[(name,strand)]={'states':int(states),'seconds':float(seconds)}
assert set(found)=={(name,strand) for name in expected for strand in '+-'}
result={}
for name in sorted(expected):
 stats=[found[(name,strand)] for strand in '+-'];result[name]={**audit['per_guide'][name],'states':sum(x['states'] for x in stats),'search_seconds':sum(x['seconds'] for x in stats),'strands':{strand:found[(name,strand)] for strand in '+-'}}
assert sum(x['states'] for x in result.values())==metrics['search_states']
print(json.dumps({'guides':result,'search_time_precision':'per-strand logged decimal values; timings include traversal and replay'},sort_keys=True))
