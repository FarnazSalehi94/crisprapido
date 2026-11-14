#!/usr/bin/env python3
"""Freeze Brunello timing inputs/commands only; never submits or runs searches.

10,000 first-in-file guides from the pinned 77,441 library; guide-order bias
explicit. Commands use one CPU, unchanged m4/b1/z1/.75/GG, independent audit,
separate both-score annotation. Stage supervisor/locked restart wrapper must
pass review before execution; the known-defective old full-run worker is not
used and is not advertised as production-ready.
"""
import argparse
import hashlib
import json
from pathlib import Path

LIBRARY_SHA256 = 'a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def select(library):
    if digest(library) != LIBRARY_SHA256: raise ValueError('library pin mismatch')
    lines = library.read_text().splitlines()
    if len(lines) != 77441*2: raise ValueError('library record count')
    pairs = list(zip(lines[::2], lines[1::2]))
    if (len({n for n,s in pairs}) != 77441 or len({s for n,s in pairs}) != 77441
            or any(not n.startswith('>') or len(s)!=20 or set(s)-set('ACGT') for n,s in pairs)):
        raise ValueError('guide identity/alphabet')
    return '\n'.join(lines[:20000])+'\n'


def prepare(root):
    meta=root/'00_manifest/timing-10000-v1'
    meta.mkdir(exist_ok=False)
    top=root.parents[1]
    library=top/'guides/brunello_77441.fa'
    guides=meta/'guides.fa'; guides.write_text(select(library))
    sbds=top/'runs/sbds_ngg_nrr_v1'
    bindings=sbds/'00_manifest/bindings/HG00735_pat'
    binding=json.loads((bindings/'index.crispr-binding.json').read_text())
    runtime=Path('/export/local/home/salehi/crisprapido-research-workflow-v1')
    repo=Path(__file__).resolve().parents[3]
    scorer=sbds/'00_manifest/nrr-tooling-v1/scorer-v2'
    out=root/'01_pilot/timing-10000-v1/HG00735_pat'
    out.mkdir(parents=True,exist_ok=False)
    cli=runtime/'target-final/release/crisprapido'; native=runtime/'pam-search'
    ref=bindings/'reference.fa'; index=bindings/'index'; ri=bindings/'reference-index.json'
    audit=repo/'scripts/release_validation/audit_emitted.py'
    samtools=Path('/home/salehi/miniconda3/envs/assembly_eval/bin/samtools')
    commands=dict(
        search=list(map(str,[cli,'--search-backend','pam-anchored','--pam-search-bin',native,
             '--columba-index',index,'-r',ref,'-G',guides,'-m',4,'-b',1,'-z',1,'-f',.75,'-p','GG','-t',1,
             '--mismatch-scores',repo/'mismatch_scores.txt','--pam-scores',repo/'pam_scores.txt'])),
        audit=list(map(str,['/usr/bin/python3',audit,'--paf',out/'eligible.cfd.paf',
             '--fasta',ref,'--guides',guides,'--matrices',repo,'--samtools',samtools,
             '--output',out/'audit','--max-mismatches',4,'--max-bulges',1,'--max-bulge-size',1])),
        score=list(map(str,[cli,'-r',ref,'-G',guides,'--score-paf',out/'eligible.cfd.paf',
             '--scoring','both','--scoring-output-dir',out/'neural',
             '--crispr-bulge-python',runtime/'inference-env/bin/python',
             '--crispr-bulge-model-root',runtime/'models-verified',
             '--scoring-script',scorer/'score_paf.py',
             '--scoring-reference-sha256',binding['reference_sha256'],
             '--scoring-reference-index',ri,'--scoring-reference-index-sha256',digest(ri),
             '--scoring-worker-timeout',43200])))
    files=[library,guides,cli,native,ri,bindings/'index.crispr-binding.json',
           repo/'mismatch_scores.txt',repo/'pam_scores.txt',audit,samtools]
    pins={str(p):digest(p) for p in files}
    pins[str(ref)]=binding['reference_sha256']
    for suffix,h in binding['components_sha256'].items(): pins[str(index)+'.'+suffix]=h
    for f in scorer.glob('*.py'): pins[str(f)]=digest(f)
    model_manifest=json.loads((scorer/'model_manifest.json').read_text())
    for model in model_manifest['models']:
        pins[str(runtime/'models-verified'/model['path'])]=model['sha256']
    plan=dict(status='PREPARED_NOT_SUBMITTED',guide_count=10000,assembly='HG00735_pat',
        configuration=dict(m=4,b=1,z=1,f=.75,p='GG'),full_cohort_authorized=False,
        commands=commands,pins=pins,output=str(out),
        resources=dict(cpus=1,mem_gib=56,time_hours=72,rss_stop_gib=48,
                       scratch_free_gb=60,output_reserve_gb=150,concurrency=1),
        estimated_seconds=dict(traversal=83158.5,audit=5360,neural=39630),
        selection='first 10000 frozen records; deterministic nonrandom guide mix',
        launch_gates=['SBDS/rebuild jobs have priority; do not compete for current resources',
                     'ensure FASTA .fai exists and matches reference; audit must not generate beside immutable reference',
                     'test supervisor/locking/checksummed stage-restart wrapper before launch',
                     'reserve 150 GB for uncompressed output + copies, keep every eligible site for IMPG'])
    (meta/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    return plan


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();plan=prepare(a.root)
    print(json.dumps(dict(status=plan['status'],guides=plan['guide_count'],output=plan['output'])))
