#!/usr/bin/env python3
"""Summarize non-overlapping 20-mer duplication; a proxy, not repeat annotation."""
import argparse
from collections import Counter
import json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--panel',required=True);p.add_argument('--output',required=True);a=p.parse_args()
root=Path(a.panel);windows=json.loads((root/'windows.json').read_text());rows=[]
for w in windows:
 seq=''.join((root/w['name']/'ref.fa').read_text().splitlines()[1:]).upper();words=[seq[i:i+20] for i in range(0,len(seq)-19,20) if set(seq[i:i+20])<=set('ACGT')]
 counts=Counter(words);rows.append({'window':w['name'],'sampled_20mers':len(words),'duplicate_occurrences':sum(n for n in counts.values() if n>1),
  'duplicate_fraction':sum(n for n in counts.values() if n>1)/len(words) if words else None,'max_frequency':max(counts.values()) if counts else None})
Path(a.output).write_text(json.dumps(rows,indent=2)+'\n')
