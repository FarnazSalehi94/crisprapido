#!/usr/bin/env python3
"""Freeze small pilot tooling/configs. Does not submit jobs or authorize full runs."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from inventory import LIBRARY_SHA256, digest, library


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', required=True, choices=['sbds', 'brunello'])
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--guides', required=True, type=Path)
    a = p.parse_args()
    a.root = a.root.resolve()
    if a.root.name != a.project+'_pam_aware_v1':
        raise ValueError('project/root mismatch')
    rows = library(a.guides)
    if len(rows) != (1 if a.project == 'sbds' else 3):
        raise ValueError('not the safety-pilot guide set')
    repo = Path(__file__).resolve().parents[2]
    runtime = Path('/export/local/home/salehi/crisprapido-research-workflow-v1')
    bound = Path('/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1')
    meta = a.root/'00_manifest'
    tooling = meta/'pilot-tooling-v1'
    tooling.mkdir()
    for name in ['worker.py', 'qc.py']:
        shutil.copyfile(Path(__file__).parent/name, tooling/name)
    shutil.copyfile(repo/'research/pam_anchored_poc/monitored_pilot/supervise.py', tooling/'supervise.py')
    guides = meta/('sbds.fa' if a.project == 'sbds' else 'brunello-pilot.fa')
    with guides.open('xb') as f:
        f.write(a.guides.read_bytes())
    plan = dict(format='hprcv2-guarded-pilot-v1', project=a.project, root=str(a.root),
                assembly='HG00735_pat', full_brunello_authorized=False, guide_count=len(rows),
                guides=str(guides), source=str(runtime/'source'),
                cli=str(runtime/'target-final/release/crisprapido'), native=str(runtime/'pam-search'),
                reference=str(bound/'reference.fa'), index=str(bound/'index'),
                reference_sha256='bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b',
                fai='/moosefs/pangenomes/HPRCv2/HG00735_pat_hprc_r2_v1.0.1.fa.gz.fai',
                reference_index='/export/local/home/salehi/crisprapido-scoring-cli-v1/reference-index.json',
                models=str(runtime/'models-verified'), model_python=str(runtime/'inference-env/bin/python'),
                samtools='/home/salehi/miniconda3/envs/assembly_eval/bin/samtools',
                supervisor=str(tooling/'supervise.py'), worker=str(tooling/'worker.py'),
                fraction=.75, pam='GG',
                configurations={'primary': [5, 1, 2], 'historical-domain-control': [2, 1, 1]} if a.project == 'sbds' else {'A': [4, 1, 1], 'B': [5, 1, 2]},
                limits=dict(rss_kib=50331648, min_start_available_kib=67108864,
                    min_running_available_kib=8388608, min_scratch_bytes=60000000000,
                    min_running_disk_bytes=5000000000, stage_timeout_seconds=6600, file_limit_bytes=10737418240),
                provenance=dict(git_commit=subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], universal_newlines=True).strip(),
                    baseline_cli_source_commit='6fa96d3', current_search_and_scoring_core='unchanged from attested development build',
                    model='Official c_2 five-classifier mean; pinned HDF5, no retraining/download/pickle inference',
                    scope='One verified assembly; three sequence-GC-stratified Brunello guides. Not a representative cohort/library estimate.'))
    pin_paths = [Path(plan[k]) for k in ['cli', 'native', 'fai', 'reference_index', 'model_python', 'samtools', 'supervisor', 'worker', 'guides']]
    pin_paths += [bound/'index.crispr-binding.json', tooling/'qc.py', runtime/'inference-env/pyvenv.cfg',
                  runtime/'source/mismatch_scores.txt', runtime/'source/pam_scores.txt',
                  runtime/'source/research/crispr_bulge/model_manifest.json']
    pin_paths += list((runtime/'source/research').rglob('*.py'))
    pin_paths += list((runtime/'source/scripts/release_validation').glob('*.py'))
    if a.project == 'brunello':
        full = Path('/moosefs/salehi/HPRCv2_CRISPRapido/guides/brunello_77441.fa')
        if digest(full) != LIBRARY_SHA256 or not set(rows) <= set(library(full)):
            raise ValueError('pilot does not exactly preserve frozen Brunello identifiers/sequences')
        plan['full_library'] = str(full)
        pin_paths.append(full)
    plan['pins'] = {str(path): digest(path) for path in sorted(set(pin_paths))}
    manifest = json.loads((runtime/'source/research/crispr_bulge/model_manifest.json').read_text())
    for model in manifest['models']:
        plan['pins'][str(runtime/'models-verified'/model['path'])] = model['sha256']
    with (meta/'pilot-plan.json').open('x') as f:
        json.dump(plan, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(plan=str(meta/'pilot-plan.json'), sha256=digest(meta/'pilot-plan.json'), worker=str(tooling/'worker.py'))))


if __name__ == '__main__':
    main()
