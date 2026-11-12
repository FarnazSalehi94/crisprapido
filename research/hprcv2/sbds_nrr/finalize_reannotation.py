#!/usr/bin/env python3
"""Publish checksummed completion only after the scoring supervisor finishes."""
import hashlib
import json
from pathlib import Path
import sys


def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def finalize(root):
    summary=json.loads((root/'summary.json').read_text())
    monitor=json.loads((root/'monitor.json').read_text())
    if summary['status']!='PASS' or monitor['reason']!='completed' or monitor['child_returncode']!=0:
        raise ValueError('incomplete scoring/validation')
    excluded=[]
    for name,sha in summary['outputs'].items():
        p=Path(name)
        # v1 included its still-open scheduler log in the inner-stage manifest.
        # Its finalized digest belongs here, after supervisor completion.
        if p.parent==root and p.name.startswith('slurm-') and p.suffix=='.log':
            excluded.append(name);continue
        if digest(p)!=sha:raise ValueError('changed output: '+name)
    files={str(p):digest(p) for p in root.rglob('*') if p.is_file() and p.name!='COMPLETE.json'}
    with (root/'COMPLETE.json').open('x') as f:
        json.dump(dict(status='PASS',outputs=files,stage_manifest_live_logs_excluded=excluded,
                       note='Finalized log hashes are included here; scientific stage artifacts all verified unchanged.'),f,indent=2)
        f.write('\n')


if __name__=='__main__':finalize(Path(sys.argv[1]).resolve())
