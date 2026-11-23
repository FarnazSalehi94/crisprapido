#!/usr/bin/env python3
"""Afterok-only publication of actual results; preserve unrelated Git changes.

No scientific work. Requires verified SUCCESS, unchanged expected branch/HEAD,
no staged changes, and the expected remote. Uses a private Git index to exclude
pre-existing handoff edits. No force push, reset, checkout, or cleanup of user data.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

HEADER = '# CRISPRapido project handoff\n'
RESULT = 'docs/BRUNELLO_10K_SCORING_RESULTS.md'
HANDOFF = 'docs/PROJECT_HANDOFF.md'


def require(ok, message):
    if not ok: raise ValueError(message)


def git(repo, *args, data=None, env=None):
    return subprocess.check_output(['git','-C',str(repo),*args],input=data,text=True,env=env).strip()


def commit_results(repo, expected_head, expected_branch, expected_remote, body, handoff_section):
    """Commit only a new result document and our prepended handoff section.

    Unstaged handoff edits remain unstaged; unrelated working files are untouched.
    Ref update is compare-and-swap. Remote push is separately fast-forward only.
    """
    require(git(repo,'branch','--show-current') == expected_branch, 'branch changed; publication stopped')
    require(git(repo,'rev-parse','HEAD') == expected_head, 'HEAD changed; publication stopped')
    require(git(repo,'remote','get-url','origin') == expected_remote, 'remote changed')
    require(not git(repo,'diff','--cached','--name-only'), 'staged changes exist; publication stopped')
    result = repo/RESULT; handoff = repo/HANDOFF
    require(not result.exists(), 'result document already exists; do not overwrite')
    before = handoff.read_text(); base = subprocess.check_output(['git','-C',str(repo),'show','HEAD:'+HANDOFF],text=True)
    require(before.startswith(HEADER) and base.startswith(HEADER), 'handoff header drift')
    working = HEADER+handoff_section+before[len(HEADER):]
    index_text = HEADER+handoff_section+base[len(HEADER):]
    # All checks above precede any working-tree modification.
    with result.open('x') as f: f.write(body)
    require(handoff.read_text() == before, 'handoff changed concurrently')
    handoff.write_text(working)
    with tempfile.TemporaryDirectory(prefix='brunello-result-index-',dir=repo/'.git') as tmp:
        env = dict(os.environ,GIT_INDEX_FILE=str(Path(tmp)/'index'))
        git(repo,'read-tree',expected_head,env=env)
        for path, text in [(RESULT,body),(HANDOFF,index_text)]:
            blob = git(repo,'hash-object','-w','--stdin',data=text,env=env)
            git(repo,'update-index','--add','--cacheinfo','100644,'+blob+','+path,env=env)
        tree = git(repo,'write-tree',env=env)
        commit = git(repo,'commit-tree',tree,'-p',expected_head,data='Record verified Brunello 10k scoring completion\n',env=env)
        git(repo,'update-ref','refs/heads/'+expected_branch,commit,expected_head)
        # Updating index entries only: never overwrite working files or stage
        # unrelated edits. These two entries now match the newly committed tree.
        for path in [RESULT,HANDOFF]:
            entry = git(repo,'ls-tree',commit,'--',path).split()
            git(repo,'update-index','--add','--cacheinfo',entry[0]+','+entry[2]+','+path)
    return commit


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--plan',required=True,type=Path); a=p.parse_args()
    plan=json.loads(a.plan.read_text()); root=Path(plan['root']); repo=Path(plan['repository'])
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == plan['publisher_sha256'], 'publisher pin drift')
    success_path=root/'SUCCESS.json'; success=json.loads(success_path.read_text())
    require(success['status']=='COMPLETE' and success['validation']['status']=='PASS' and success['validation']['rows']==57578781, 'scientific completion not verified')
    require(success['scoring_job']==plan['scoring_job'] and success['plan_sha256']==plan['scoring_plan_sha256'], 'stage identity drift')
    require(success['supervisor']['reason']=='completed' and success['supervisor']['child_returncode']==0, 'supervisor not successful')
    require((root/'RESULTS.md').is_file(), 'finalizer results absent')
    digest=hashlib.sha256(success_path.read_bytes()).hexdigest()
    v=success['validation']; memory=success['supervisor']; availability=json.dumps(v['availability'],sort_keys=True)
    body=(f'# Brunello 10k scoring verified complete\n\n'
          f'Scoring job **{success["scoring_job"]} COMPLETED 0:0**; afterok finalizer **{success["finalizer_job"]}**.\n\n'
          f'Independent full verification PASS: **{v["rows"]:,} rows**, exact source coverage/order/CFD and all non-CB annotations, '
          f'CB availability unchanged from the hash-verified historical conversion spool: `{availability}`. '
          'Historical predictions were not reused.\n\n'
          f'Supervisor elapsed {memory["elapsed_seconds"]:.1f} seconds; sampled summed-tree RSS peak '
          f'{memory["sampled_peak_rss_kib"]:,} KiB (not unique physical memory; sampling may miss spikes and double-count shared pages).\n\n'
          f'Full evidence: `{root}`. SUCCESS SHA256 `{digest}`. Frozen scoring plan SHA256 `{success["plan_sha256"]}`.\n\n'
          'The pre-run 101k comparison max CB delta was 1.51e-9 within predeclared 1e-6 tolerance. '
          'Global duplicate checking and content-bound resumable chunking precede verified merge. Original attempts and search/audit checkpoints remain intact. '
          '**No SBDS, new search, full library or full cohort run.**\n')
    section=(f'\n## {datetime.now(timezone.utc).date()} Brunello 10k scoring completion — verified after long run\n\n'
             f'Scoring **{success["scoring_job"]}**, afterok finalizer **{success["finalizer_job"]}**: '
             f'new-stage SUCCESS / full-stream validation PASS, {v["rows"]:,} rows. '
             f'Actual results: [BRUNELLO_10K_SCORING_RESULTS.md](BRUNELLO_10K_SCORING_RESULTS.md). '
             f'Evidence `{root}`. All prior failures/checkpoints preserved; no SBDS/search/cohort run.\n')
    commit=commit_results(repo,plan['expected_head'],plan['branch'],plan['remote'],body,section)
    # No force push: concurrent remote advancement must fail visibly.
    subprocess.run(['git','-C',str(repo),'push','origin',plan['branch']],check=True)
    (root/'PUBLICATION.json').write_text(json.dumps(dict(status='PUSHED',commit=commit,branch=plan['branch'],success_sha256=digest),indent=2)+'\n')
    (root/'REPORT.md').write_text(body+'\n---\n\n'+(root/'REPORT.md').read_text())
    (root/'PROJECT_HANDOFF.md').write_text((repo/HANDOFF).read_text())
    print(json.dumps(dict(status='PUSHED',commit=commit)))


if __name__=='__main__':
    main()
