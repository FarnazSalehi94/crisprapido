#!/usr/bin/env python3
"""Validate retained metadata evidence; does not open HPRC indexes."""
import json
from pathlib import Path
p=Path(__file__).with_name('results');read=lambda x:json.loads((p/x).read_text());panel=read('panel_metadata.json');aliases=read('alias_metadata.json');resolved=read('resolved_metadata.json');summary=read('summary.json');gate=read('direct_gate.json')
assert len(panel)==summary['panel_rows']==464 and len({x['assembly'] for x in panel})==464
ok=[x for x in panel if x['status']=='metadata_observed'];assert len(ok)==459 and all(x['metadata']=='21\n4\nVANILLA\n' and len(x['sizes'])==12 and all(x['sizes'].values()) and not x['binding_exists'] for x in ok)
assert len(aliases)==4 and all(x['status']=='metadata_observed' and x['meta']=='21|4|VANILLA' and x['dictionary_matches'] and x['offsets_match'] and not x['binding_exists'] for x in aliases)
assert len(resolved)==12 and all(x['status']=='metadata_observed' and x['dictionary_matches'] and x['offsets_match'] and x['text_header_matches'] and not x['binding_exists'] for x in resolved)
assert gate['returncode']!=0 and gate['stdout_bytes']==0 and 'manifest missing' in gate['stderr']
print('HPRC metadata inventory, aliases, resolved set and fail-closed gate consistent')
