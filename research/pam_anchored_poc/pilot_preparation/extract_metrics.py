#!/usr/bin/env python3
"""Require timed-native metrics and report audited output counts."""
import json,re,sys
from pathlib import Path
log=Path(sys.argv[1]).read_text();audit=json.loads(Path(sys.argv[2]).read_text());idx=re.findall(r'^index_load_seconds=([0-9.eE+-]+)$',log,re.M);search=re.findall(r'^states=(\d+) search_wall_seconds=([0-9.eE+-]+)$',log,re.M)
assert len(idx)==len(search)==1 and audit['status']=='pass' and audit['known_exact_target_recovered']
print(json.dumps({'native_index_load_seconds':float(idx[0]),'native_search_states':int(search[0][0]),'native_search_wall_seconds':float(search[0][1]),'reported_sites':audit['reported'],'known_exact_target_recovered':True},sort_keys=True))
