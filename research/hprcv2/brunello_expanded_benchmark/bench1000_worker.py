#!/usr/bin/env python3
"""Bounded 1,000-guide scaling benchmark at the recommended configuration.

Same guards/stages as the 16-guide benchmark (search -> audit -> both-score
annotation -> summary), one invocation per assembly. The panel is the first
1,000 records of the frozen library in file order: deterministic, not random.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
for candidate in (HERE, HERE.parent, HERE.parent.parent):
    sys.path.insert(0, str(candidate))
from bench16_worker import (CONFIGS, LIMITS, digest, inner as pilot_inner,  # noqa: E402
                            require, save, validate)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', required=True, type=Path)
    p.add_argument('--assembly', required=True)
    p.add_argument('--config', required=True)
    p.add_argument('--inner', type=Path)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    require(plan['format'] == 'brunello-1000-benchmark-v1', 'plan format drift')
    require(plan['guide_count'] == 1000 and plan['full_run_authorized'] is False, 'scope drift')
    require(a.config in CONFIGS and CONFIGS[a.config] == plan['config_parameters'], 'configuration drift')
    require(a.assembly in plan['assemblies'], 'unknown assembly')
    for assembly, spec in plan['assemblies'].items():
        binding = json.loads(Path(spec['binding']).read_text())
        require(binding['reference_sha256'] == spec['reference_sha256'], 'binding mismatch: '+assembly)
        require(Path(spec['scoring_reference_index']).is_absolute()
                and spec['scoring_reference_index'] in plan['pins'], 'unpinned scoring index')
    plan['_plan_sha256'] = digest(a.plan)
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURMD_NODENAME'), 'Slurm execution required')
    require(int(os.environ['SLURM_CPUS_PER_TASK']) >= 8, 'reserve at least eight CPUs; application uses one')
    if a.inner:
        pilot_inner(dict(plan, guide_count=1000, limits=LIMITS, configs=CONFIGS), a.assembly, a.config, a.inner)
        return
    logical = Path(plan['root'])/'01_pilot'/'benchmark-1000-v1'/(a.assembly+'-'+a.config)
    logical.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (logical/'task.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        success = logical/'SUCCESS.json'
        if success.exists():
            result = json.loads(success.read_text())
            require(result['plan_sha256'] == plan['_plan_sha256'], 'completed task has a different plan')
            for path, spec in result['outputs'].items():
                require(digest(path) == spec['sha256'], 'completed output changed: '+path)
            print('Already complete; checksums verified; no rerun.')
            return
        attempt = logical/('attempt-'+os.environ['SLURM_JOB_ID'])
        attempt.mkdir()
        save(attempt/'INCOMPLETE.json', dict(plan=str(a.plan), assembly=a.assembly, config=a.config))
        try:
            require(all(digest(path) == expected for path, expected in plan['pins'].items()), 'pin mismatch')
            from bench16_worker import verify_pins
            verify_pins(plan)
            import subprocess
            command = list(map(str, ['/usr/bin/python3', plan['supervisor'], '--limit-kib', LIMITS['rss_kib'],
                       '--interval', '.1', '--log', attempt/'monitor.log', '--summary', attempt/'monitor-summary.json',
                       '--', '/usr/bin/python3', Path(__file__).resolve(), '--plan', a.plan,
                       '--assembly', a.assembly, '--config', a.config, '--inner', attempt]))
            subprocess.run(command, check=True)
            monitor = json.loads((attempt/'monitor-summary.json').read_text())
            require(monitor['reason'] == 'completed', 'supervisor did not pass')
            data = json.loads((attempt/'STAGES_COMPLETE.json').read_text())
            for path, spec in data['outputs'].items():
                require(digest(path) == spec['sha256'], 'output checksum mismatch')
            for name in ['monitor.log', 'monitor-summary.json', 'STAGES_COMPLETE.json']:
                path = attempt/name
                data['outputs'][str(path)] = dict(sha256=digest(path), bytes=path.stat().st_size)
            save(success, dict(status='COMPLETE', plan_sha256=plan['_plan_sha256'], attempt=str(attempt),
                              outputs=data['outputs'], monitor=monitor))
            (attempt/'INCOMPLETE.json').unlink()
        except BaseException as exc:
            save(attempt/'FAILED.json', dict(error=str(exc), plan_sha256=plan['_plan_sha256']))
            raise


if __name__ == '__main__':
    main()
