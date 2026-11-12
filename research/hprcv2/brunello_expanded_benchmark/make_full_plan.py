#!/usr/bin/env python3
"""Build the full-run plan AFTER explicit user authorization. Launches nothing.

Usage (only after the user approves the configuration and scope):

    python3 make_full_plan.py --root ... --bindings registry.json \
        --config candidate-A --approved-by '<user decision>' --decision-ref '<reference>'

The authorization record is copied verbatim into the plan; this script cannot
invent it. SBDS is paused and is not part of this plan in any way.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from inventory import LIBRARY_SHA256, digest

CONFIGS = dict(historical=[2, 1, 1], intermediate=[3, 1, 1],
               **{'candidate-A': [4, 1, 1], 'candidate-B': [5, 1, 2]})
SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
LIMITS = dict(chunk_timeout_seconds=86400, scoring_timeout_seconds=86400,
              rss_kib=75161927, min_running_available_kib=8388608,
              min_running_scratch_bytes=100000000000, min_running_root_bytes=200000000000,
              file_limit_bytes=107374182400)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--bindings', required=True, type=Path, help='registry JSON from prepare_binding/import steps')
    p.add_argument('--config', required=True, choices=CONFIGS)
    p.add_argument('--chunk-size', type=int, required=True)
    p.add_argument('--approved-by', required=True)
    p.add_argument('--decision-ref', required=True)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    if a.root.name != 'brunello_pam_aware_v1':
        raise ValueError('unexpected project root')
    if not (1000 <= a.chunk_size <= 77441):
        raise ValueError('chunk size must bound restart loss while amortizing initialization')
    registry = json.loads(a.bindings.read_text())
    repo = Path(__file__).resolve().parents[3]
    runtime = Path('/export/local/home/salehi/crisprapido-research-workflow-v1')
    meta = a.root/'00_manifest'
    tooling = meta/'full-run-tooling-v1'
    tooling.mkdir(exist_ok=True)
    for name in ['full_run_worker.py', 'verify_full_run.py', 'submit_full_run.slurm']:
        shutil.copyfile(Path(__file__).parent/name, tooling/name)
    for name in ['qc.py', 'collect.py', 'inventory.py']:
        shutil.copyfile(Path(__file__).resolve().parents[1]/name, tooling/name)
    shutil.copyfile(repo/'research/pam_anchored_poc/monitored_pilot/supervise.py', tooling/'supervise.py')
    guides = Path('/moosefs/salehi/HPRCv2_CRISPRapido/guides/brunello_77441.fa')
    plan = dict(format='brunello-full-run-v1', root=str(a.root), guides=str(guides),
                guide_count=77441, config=CONFIGS[a.config], chunk_size=a.chunk_size,
                authorization=dict(approved_by=a.approved_by, decision_ref=a.decision_ref,
                                   scope='Brunello 77,441 guides x 464 HPRCv2 haplotype assemblies; SBDS excluded (paused pending PAM confirmation)'),
                limits=LIMITS,
                cli=str(runtime/'target-final/release/crisprapido'), native=str(runtime/'pam-search'),
                matrices=dict(mismatch=str(repo/'mismatch_scores.txt'), pam=str(repo/'pam_scores.txt')),
                release_scripts=str(repo/'scripts/release_validation'),
                samtools='/home/salehi/miniconda3/envs/assembly_eval/bin/samtools',
                models=str(runtime/'models-verified'), model_python=str(runtime/'inference-env/bin/python'),
                scoring_script=str(runtime/'source/research/crispr_bulge/score_paf.py'),
                supervisor=str(tooling/'supervise.py'), worker=str(tooling/'full_run_worker.py'),
                assemblies=registry['assemblies'],
                provenance=dict(git_commit=subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], universal_newlines=True).strip(),
                                library_sha256=LIBRARY_SHA256,
                                note='Every assembly binding must come from prepare_binding.py or the trusted import path; no metadata-only acceptance.'))
    pins = [Path(plan[k]) for k in ['cli', 'native', 'samtools', 'model_python', 'scoring_script', 'supervisor', 'worker', 'guides']]
    pins += [Path(plan['matrices'][k]) for k in ['mismatch', 'pam']]
    pins += [Path(plan['release_scripts'])/'audit_emitted.py']
    for spec in plan['assemblies'].values():
        pins += [Path(spec[k]) for k in ['reference', 'binding', 'fai', 'scoring_reference_index']]
        pins += [Path(str(spec['index'])+'.'+s) for s in SUFFIXES]
    plan['pins'] = {str(p): digest(p) for p in sorted(set(pins))}
    manifest = json.loads((runtime/'source/research/crispr_bulge/model_manifest.json').read_text())
    for model in manifest['models']:
        plan['pins'][str(runtime/'models-verified'/model['path'])] = model['sha256']
    with a.out.open('x') as f:
        json.dump(plan, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(plan=str(a.out), sha256=digest(a.out), assemblies=len(plan['assemblies']),
                          config=a.config, chunk_size=a.chunk_size)))


if __name__ == '__main__':
    main()
