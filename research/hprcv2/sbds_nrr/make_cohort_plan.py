#!/usr/bin/env python3
"""Freeze a restartable SBDS NGG/NRR cohort plan from completed bindings.

Includes ONLY assemblies with a committed binding.json. Excluded assemblies
are recorded explicitly; the plan never guesses or zero-fills them. Requires
an explicit authorization token recorded verbatim from the user for any scope
beyond the controlled test assemblies.
"""
import argparse
import hashlib
import json
from pathlib import Path

SHARED_FILES = ['scripts/release_validation/audit_emitted.py']


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--repo', type=Path, required=True)
    p.add_argument('--bindings', type=Path, required=True)
    p.add_argument('--include', required=True, help='comma-separated assemblies to include')
    p.add_argument('--authorization', required=True, help='verbatim user authorization for this scope')
    p.add_argument('--scorer-dir', type=Path, required=True, help='frozen NRR-enabled scorer bundle directory')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    source = Path('/export/local/home/salehi/crisprapido-research-workflow-v1/source')
    runtime = Path('/export/local/home/salehi/crisprapido-research-workflow-v1')
    pilot = json.loads(Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1/00_manifest/pilot-plan-v1.json').read_text())
    include = a.include.split(',')
    assemblies = {}
    for name in include:
        stage = a.bindings/name
        binding = json.loads((stage/'binding.json').read_text())
        if (stage/'INCOMPLETE').exists():
            raise ValueError('binding not committed: '+name)
        index = str(stage/'index')
        reference = str(stage/'reference.fa')
        ri = stage/'reference-index.json'
        assemblies[name] = dict(assembly=name, cli=pilot['cli'], native=pilot['native'],
            index=index, reference=reference, reference_sha256=binding['reference_sha256'],
            fai=reference+'.fai', reference_index=str(ri),
            pins={reference: binding['reference_sha256'],
                  reference+'.fai': digest(reference+'.fai'),
                  index+'.crispr-binding.json': digest(index+'.crispr-binding.json'),
                  str(ri): digest(ri)})
    shared = {str(a.repo/f): digest(a.repo/f) for f in SHARED_FILES}
    shared[str(pilot['guides'])] = digest(pilot['guides'])
    shared[str(a.bindings.parent/'nrr-tooling-v1/cohort_worker.py')] = digest(a.bindings.parent/'nrr-tooling-v1/cohort_worker.py')
    shared[str(a.bindings.parent/'nrr-tooling-v1/supervise.py')] = digest(a.bindings.parent/'nrr-tooling-v1/supervise.py')
    shared[str(source/'mismatch_scores.txt')] = digest(source/'mismatch_scores.txt')
    shared[str(source/'pam_scores.txt')] = digest(source/'pam_scores.txt')
    shared[str(source/'research/crispr_bulge/score_paf.py')] = digest(source/'research/crispr_bulge/score_paf.py')
    for f in sorted(a.scorer_dir.glob('*.py'))+[a.scorer_dir/'model_manifest.json']:
        shared[str(f)] = digest(f)
    plan = dict(format='sbds-nrr-cohort-v1',
        scoring_script=str(a.scorer_dir/'score_paf.py'),
        root='/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1',
        source=str(source), guides=pilot['guides'], models=pilot['models'],
        model_python=pilot['model_python'], samtools=pilot['samtools'],
        auditor=str(a.repo/'scripts/release_validation/audit_emitted.py'),
        assemblies=assemblies, pins=shared,
        authorization=a.authorization,
        limits=dict(rss_kib=50331648, min_start_available_kib=67108864,
            min_running_available_kib=8388608, min_scratch_bytes=60000000000,
            min_running_disk_bytes=5000000000, stage_timeout_seconds=14400, file_limit_bytes=10737418240),
        excluded=dict(note='assemblies without committed bindings are never searched or zero-filled'),
        scope='One-guide four-PAM per assembly; search m5/b1/z2/f0.75; CFD and NRR-enabled CRISPR-Bulge with explicit N/A; noncanonical biological validity not established.')
    a.out.write_text(json.dumps(plan, indent=2)+'\n')
    print(json.dumps(dict(plan=str(a.out), plan_sha256=digest(a.out), assemblies=sorted(assemblies))))


if __name__ == '__main__':
    main()
