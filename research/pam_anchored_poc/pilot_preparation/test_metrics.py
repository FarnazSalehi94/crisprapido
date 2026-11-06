#!/usr/bin/env python3
import json,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as d:
 d=Path(d);log=d/'log';audit=d/'audit';audit.write_text(json.dumps({'status':'pass','known_exact_target_recovered':True,'reported':172}));cmd=[sys.executable,str(root/'extract_metrics.py'),str(log),str(audit)]
 log.write_text('index_load_seconds=0.12\nstates=42 search_wall_seconds=1.5\n');r=subprocess.run(cmd,check=True,capture_output=True,text=True);data=json.loads(r.stdout);assert data['native_search_states']==42 and data['reported_sites']==172
 log.write_text('states=42 search_wall_seconds=1.5\n');assert subprocess.run(cmd,capture_output=True).returncode!=0
print('metric extraction passes complete logs and rejects missing index timing')
