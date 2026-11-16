#!/usr/bin/env python3
"""Freeze one disjoint recovery wave and guarded wrappers. Never submits jobs.

Selection uses an explicitly supplied manifest, excludes every committed binding
and existing staged/promoted set, and refuses overlap with supplied prior lists.
The parent must check live scheduler state and approve resource/dependency flags.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def choose(rows, root, excluded, count):
    names = [r['assembly'] for r in rows]
    if len(names) != 466 or len(set(names)) != 466:
        raise ValueError('expected exactly 466 unique assemblies')
    selected = []
    for row in rows:
        name = row['assembly']
        if '/' in name or name in ('', '.', '..'):
            raise ValueError('invalid assembly')
        if row['status'] not in ('INDEX_EIO', 'INDEX_FAILED_MISSING_INDEX', 'MISSING_INDEX'):
            continue
        binding = root/'runs/sbds_ngg_nrr_v1/00_manifest/bindings'/name
        if (name in excluded or (binding/'binding.json').exists()
                or (root/'indices_rebuilt_v1'/name).exists()
                or (root/'indices'/name/'promotion.json').exists()):
            continue
        selected.append(name)
    if len(selected) < count:
        raise ValueError('fewer eligible assemblies than requested; review smaller wave')
    return selected[:count]


def load_pins(template):
    record = template/'tooling-sha256.json'
    if record.is_file():
        return {str(template/name): sha for name, sha in json.loads(record.read_text()).items()}
    checksums = template/'tooling.sha256'
    if not checksums.is_file():
        raise ValueError('template has neither tooling-sha256.json nor tooling.sha256')
    pins = {}
    for line in checksums.read_text().splitlines():
        sha, name = line.split(None, 1)
        path = template/name.strip()
        if not path.is_file():
            raise ValueError('pinned template file missing: '+name)
        if digest(path) != sha:
            raise ValueError('template pin mismatch: '+name)
        pins[str(path)] = sha
    return pins


def prepare(root, manifest_path, template, dest, count, exclude):
    if not 1 <= count <= 48:
        raise ValueError('controlled wave size must be 1..48')
    manifest = json.loads(manifest_path.read_text())
    if manifest['cohort_size'] != 466 or not manifest['mutually_exclusive']:
        raise ValueError('invalid cohort manifest')
    excluded = set()
    for path in exclude:
        excluded.update(path.read_text().splitlines())
    names = choose(manifest['rows'], root, excluded, count)
    pins = load_pins(template)
    for path, sha in pins.items():
        if Path(path).parent != template or digest(path) != sha:
            raise ValueError('template pin/path mismatch: ' + path)
    dest.mkdir(exist_ok=False)
    for path in pins:
        source = Path(path)
        content = source.read_text().replace(str(template), str(dest)).replace('wave2-', dest.name+'-')
        (dest/source.name).write_text(content)
    (dest/'wave_guard.py').write_text(Path(__file__).with_name('wave_guard.py').read_text())
    (dest/'assemblies.txt').write_text(''.join(n+'\n' for n in names))
    for stage, limit in [('build', 60817408), ('recheck', 6291456), ('promote', 7340032)]:
        text = f'''#!/bin/bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
D={dest}
export ASSEMBLIES_FILE="$D/assemblies.txt"
exec python3 "$D/wave_guard.py" --wave-dir "$D" -- python3 "$D/supervise.py" --limit-kib {limit} --interval .5 --log "$D/{stage}-${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}.jsonl" --summary "$D/{stage}-${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}.json" -- bash "$D/{stage}-checked.sh"
'''
        (dest/(stage+'-guarded.slurm')).write_text(text)
        (dest/(stage+'-checked.sh')).write_text(f'''#!/bin/bash
set -euo pipefail
cd {dest}
# Run inside the guard so an integrity failure latches shared STOP.
sha256sum --check tooling.sha256
exec bash {stage}.slurm
''')
    frozen = {p.name: digest(p) for p in sorted(dest.iterdir()) if p.is_file()}
    (dest/'tooling.sha256').write_text(''.join(f'{sha}  {name}\n' for name, sha in frozen.items()))
    record = dict(manifest=str(manifest_path), manifest_sha256=digest(manifest_path),
        template=str(template), template_pins=pins, assemblies=names,
        excluded_lists={str(p): digest(p) for p in exclude}, frozen=frozen,
        scope='PREPARED_ONLY; main agent must review live queue/resources and submit',
        stop_policy='any guarded stage failure latches STOP; no automatic retries/deletions')
    (dest/'preparation.json').write_text(json.dumps(record, indent=2)+'\n')
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--template', type=Path, required=True)
    p.add_argument('--destination', type=Path, required=True)
    p.add_argument('--count', type=int, default=48)
    p.add_argument('--exclude', type=Path, action='append', default=[])
    a = p.parse_args()
    print(json.dumps(prepare(a.root, a.manifest, a.template, a.destination, a.count, a.exclude)))
