#!/usr/bin/env python3
"""Freeze the bounded 1,000-guide benchmark plan for one chosen configuration."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from inventory import digest

CONFIGS = dict(historical=[2, 1, 1], intermediate=[3, 1, 1],
               **{'candidate-A': [4, 1, 1], 'candidate-B': [5, 1, 2]})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--config', required=True, choices=CONFIGS)
    p.add_argument('--guides', required=True, type=Path)
    a = p.parse_args()
    a.root = a.root.resolve()
    if a.root.name != 'brunello_pam_aware_v1':
        raise ValueError('unexpected project root')
    if len([line for line in a.guides.read_text().splitlines() if line.startswith('>')]) != 1000:
        raise ValueError('guides file must contain exactly 1,000 records')
    repo = Path(__file__).resolve().parents[3]
    runtime = Path('/export/local/home/salehi/crisprapido-research-workflow-v1')
    meta = a.root/'00_manifest'
    tooling = meta/'expanded-benchmark-tooling-v1'
    for name in ['bench1000_worker.py', 'bench16_worker.py', 'analyze16.py']:
        shutil.copyfile(Path(__file__).parent/name, tooling/name)
    assemblies = {}
    for assembly, stage, scoring_index in [
            ('HG00735_pat', '/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1',
             '/export/local/home/salehi/crisprapido-scoring-cli-v1/reference-index.json'),
            ('HG01496_pat', '/export/local/home/salehi/crisprapido-hprc-import-hg01496-20261112-v1', None)]:
        bound = Path(stage)
        assemblies[assembly] = dict(reference=str(bound/'reference.fa'), index=str(bound/'index'),
                                    binding=str(bound/'index.crispr-binding.json'),
                                    fai='/moosefs/pangenomes/HPRCv2/%s_hprc_r2_v1.0.1.fa.gz.fai' % assembly,
                                    scoring_reference_index=scoring_index or str(bound/'reference-index.json'),
                                    reference_sha256=json.loads((bound/'import.json').read_text())['reference_sha256']
                                    if (bound/'import.json').exists() else
                                    json.loads((bound/'binding.json').read_text())['reference_sha256'])
    plan = dict(format='brunello-1000-benchmark-v1', root=str(a.root), full_run_authorized=False,
                guide_count=1000, guides=str(a.guides), config_parameters=CONFIGS[a.config],
                assemblies=assemblies, limits=dict(stage_timeout_seconds=14400, min_start_available_kib=67108864,
                    min_running_available_kib=8388608, min_scratch_bytes=60000000000,
                    min_running_disk_bytes=20000000000, file_limit_bytes=107374182400, rss_kib=50331648),
                cli=str(runtime/'target-final/release/crisprapido'), native=str(runtime/'pam-search'),
                matrices=dict(mismatch=str(repo/'mismatch_scores.txt'), pam=str(repo/'pam_scores.txt')),
                release_scripts=str(repo/'scripts/release_validation'),
                samtools='/home/salehi/miniconda3/envs/assembly_eval/bin/samtools',
                models=str(runtime/'models-verified'), model_python=str(runtime/'inference-env/bin/python'),
                scoring_script=str(runtime/'source/research/crispr_bulge/score_paf.py'),
                supervisor=str(tooling/'supervise.py'),
                provenance=dict(git_commit=subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], universal_newlines=True).strip(),
                    scope='Bounded 1,000-guide scaling benchmark (first 1,000 frozen-library records) on the two verified assemblies; not a full-library or cohort run.',
                    selection_basis='First records in frozen file order; deterministic; guide-mix bias unmeasured.'))
    pins = [Path(plan[k]) for k in ['cli', 'native', 'samtools', 'model_python', 'scoring_script', 'supervisor', 'guides']]
    pins += [Path(plan['matrices'][k]) for k in ['mismatch', 'pam']]
    pins += [Path(plan['release_scripts'])/'audit_emitted.py',
             tooling/'bench1000_worker.py', tooling/'bench16_worker.py', tooling/'qc.py']
    for spec in assemblies.values():
        pins += [Path(spec[k]) for k in ['reference', 'binding', 'fai', 'scoring_reference_index']]
        pins += [Path(str(spec['index'])+'.'+s) for s in
                 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()]
    plan['pins'] = {str(p): digest(p) for p in sorted(set(pins))}
    manifest = json.loads((runtime/'source/research/crispr_bulge/model_manifest.json').read_text())
    for model in manifest['models']:
        plan['pins'][str(runtime/'models-verified'/model['path'])] = model['sha256']
    path = meta/('benchmark1000-plan-%s.json' % a.config)
    with path.open('x') as f:
        json.dump(plan, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(plan=str(path), sha256=digest(path), config=a.config)))


if __name__ == '__main__':
    main()
