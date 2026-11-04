#!/usr/bin/env python3
"""Check enumerator against selected guides from an independent 172-site panel."""
import json,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[3];panel=root/'research/pam_anchored_poc/validation/results'
ids=['CHM13_nonrepeat_candidate_plus_exact','CHM13_nonrepeat_candidate_minus_mm5','CHM13_nonrepeat_candidate_plus_I2mm3','CHM13_repeat_candidate_minus_D2mm5','GRCh38_repeat_candidate_plus_D2','GRCh38_nonrepeat_candidate_minus_I1']
lines=(panel/'guides.fa').read_text().splitlines();guides={line[1:]:lines[i+1] for i,line in enumerate(lines) if line.startswith('>')};assert all(x in guides for x in ids)
with tempfile.TemporaryDirectory() as td:
 d=Path(td);(d/'guide.fa').write_text(''.join('>'+x+'\n'+guides[x]+'\n' for x in ids));(d/'empty.paf').write_text('')
 subprocess.run([sys.executable,str(Path(__file__).with_name('audit.py')),'--reference',str(panel/'ref.fa'),'--guides',str(d/'guide.fa'),'--baseline',str(d/'empty.paf'),'--experimental',str(d/'empty.paf'),'--out',str(d/'summary.json'),'--oracle'],check=True,stdout=subprocess.DEVNULL)
 got={tuple(x) for x in json.loads((d/'summary.json').read_text())['oracle']['baseline_missing']}
 expected={(x['guide'],x['contig'],x['strand'],x['anchor']+(3 if x['strand']=='-' else 0)) for x in map(json.loads,(panel/'oracle.jsonl').open()) if x['guide'] in ids}
 assert got==expected,(sorted(got-expected),sorted(expected-got))
 print('%d-guide existing 172-panel oracle: %d/%d anchors' %(len(ids),len(got&expected),len(expected)))
