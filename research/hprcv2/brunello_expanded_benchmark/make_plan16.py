#!/usr/bin/env python3
"""Freeze the expanded-benchmark plan; no submission, no full-library scope."""
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
LIMITS = dict(stage_timeout_seconds=10800, min_start_available_kib=67108864,
              min_running_available_kib=8388608, min_scratch_bytes=60000000000,
              min_running_disk_bytes=10000000000, file_limit_bytes=21474836480,
              rss_kib=50331648)
BOUNDS = {
    'HG00735_pat': dict(root='/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1',
                        fai='/moosefs/pangenomes/HPRCv2/HG00735_pat_hprc_r2_v1.0.1.fa.gz.fai',
                        reference_sha256='bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b',
                        scoring_reference_index='/export/local/home/salehi/crisprapido-scoring-cli-v1/reference-index.json'),
    'HG01496_pat': dict(root='/export/local/home/salehi/crisprapido-hprc-import-hg01496-20261112-v1',
                        fai='/moosefs/pangenomes/HPRCv2/HG01496_pat_hprc_r2_v1.0.1.fa.gz.fai',
                        reference_sha256='f69e5392fb0bc48807c11d439471ead7832486e092a3d67123d79d399e4a40c8',
                        scoring_reference_index='/export/local/home/salehi/crisprapido-hprc-import-hg01496-20261112-v1/reference-index.json')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--guides', required=True, type=Path)
    a = p.parse_args()
    a.root = a.root.resolve()
    if a.root.name != 'brunello_pam_aware_v1':
        raise ValueError('unexpected project root')
    repo = Path(__file__).resolve().parents[3]
    runtime = Path('/export/local/home/salehi/crisprapido-research-workflow-v1')
    meta = a.root/'00_manifest'
    tooling = meta/'expanded-benchmark-tooling-v1'
    tooling.mkdir()
    for name in ['bench16_worker.py', 'analyze16.py', 'test_expanded.py']:
        shutil.copyfile(Path(__file__).parent/name, tooling/name)
    for name in ['qc.py', 'collect.py', 'inventory.py', 'compare_sbds.py']:
        shutil.copyfile(Path(__file__).resolve().parents[1]/name, tooling/name)
    shutil.copyfile(repo/'research/pam_anchored_poc/monitored_pilot/supervise.py', tooling/'supervise.py')
    assemblies = {}
    for name, bound in BOUNDS.items():
        assemblies[name] = dict(reference=str(Path(bound['root'])/'reference.fa'),
                                 reference_sha256=bound['reference_sha256'],
                                 index=str(Path(bound['root'])/'index'),
                                 binding=str(Path(bound['root'])/'index.crispr-binding.json'),
                                 fai=bound['fai'], scoring_reference_index=bound['scoring_reference_index'])
    plan = dict(format='brunello-expanded-benchmark-v1', root=str(a.root),
                full_run_authorized=False, guide_count=16, guides=str(a.guides),
                assemblies=assemblies, configs=CONFIGS, limits=LIMITS,
                cli=str(runtime/'target-final/release/crisprapido'), native=str(runtime/'pam-search'),
                matrices=dict(mismatch=str(repo/'mismatch_scores.txt'), pam=str(repo/'pam_scores.txt')),
                release_scripts=str(repo/'scripts/release_validation'),
                samtools='/home/salehi/miniconda3/envs/assembly_eval/bin/samtools',
                models=str(runtime/'models-verified'), model_python=str(runtime/'inference-env/bin/python'),
                scoring_script=str(runtime/'source/research/crispr_bulge/score_paf.py'),
                supervisor=str(tooling/'supervise.py'), worker=str(tooling/'bench16_worker.py'),
                provenance=dict(git_commit=subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], universal_newlines=True).strip(),
                    panel='16 deterministic guides: 15 historical-retained-count x GC cells plus one non-targeting control; median-count selection, timing-blind.',
                    scope='Bounded 16-guide benchmark across two verified assemblies and four configurations; no full-library or cohort launch.',
                    f_match_note='f=0.75 requires >=15/20 matched bases; all four configurations permit their full mismatch+bulge domains without clipping.'))
    pins = [Path(plan[k]) for k in ['cli', 'native', 'samtools', 'model_python', 'scoring_script', 'supervisor', 'worker', 'guides']]
    pins += [Path(plan['matrices'][k]) for k in ['mismatch', 'pam']]
    pins += [Path(plan['release_scripts'])/'audit_emitted.py', tooling/'qc.py', tooling/'analyze16.py']
    for spec in assemblies.values():
        pins += [Path(spec[k]) for k in ['reference', 'index.crispr-binding.json' if False else 'binding', 'fai', 'scoring_reference_index']]
        pins += [Path(str(spec['index'])+'.'+s) for s in
                 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()]
    plan['pins'] = {str(p): digest(p) for p in sorted(set(pins))}
    manifest = json.loads((runtime/'source/research/crispr_bulge/model_manifest.json').read_text())
    for model in manifest['models']:
        plan['pins'][str(runtime/'models-verified'/model['path'])] = model['sha256']
    path = meta/'expanded-benchmark-plan.json'
    with path.open('x') as f:
        json.dump(plan, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(plan=str(path), sha256=digest(path))))


if __name__ == '__main__':
    main()
