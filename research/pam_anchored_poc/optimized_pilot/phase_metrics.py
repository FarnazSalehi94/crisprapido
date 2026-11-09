#!/usr/bin/env python3
"""Extract instrumentation without claiming unavailable output-only timing."""
import json,re,sys
from pathlib import Path
log=Path(sys.argv[1]).read_text();timing=dict(x.split('=',1) for x in Path(sys.argv[2]).read_text().splitlines());audit=json.loads(Path(sys.argv[3]).read_text());fields={}
for name in ('snapshot_copy_seconds','snapshot_verify_pre_seconds','native_process_seconds','native_binding_seconds','native_reference_load_seconds','native_index_load_seconds','snapshot_verify_post_seconds','python_reference_load_seconds','runner_total_seconds'):
 vals=re.findall(r'^PROFILE '+name+r'=([0-9.eE+-]+)$',log,re.M);assert len(vals)==1,(name,len(vals));fields[name]=float(vals[0])
m=re.findall(r'^states=(\d+) search_wall_seconds=([0-9.eE+-]+)$',log,re.M);assert len(m)==1;fields['search_states']=int(m[0][0]);fields['native_search_seconds_including_replay']=float(m[0][1]);fields['workflow_seconds_whole_seconds']=int(timing['elapsed_seconds']);fields['workflow_exit_code']=int(timing['exit_code']);fields['reported_sites']=audit['reported'];fields['combined_rust_startup_parse_scoring_format_output_seconds_proxy']=fields['workflow_seconds_whole_seconds']-fields['runner_total_seconds'];fields['qualification']='combined residual is not isolated CFD or output; Bash elapsed has one-second granularity; residual can be negative due to rounding'
assert fields['workflow_exit_code']==0 and audit['status']=='pass';print(json.dumps(fields,sort_keys=True))
