#!/usr/bin/env python3
"""Scoring-only migration for the existing Brunello 10k pilot; no search code.

prepare freezes a separate stage identity; work verifies and reuses the original
search/audit, then runs sequential bounded chunks. finalize requires an afterok
Slurm dependency and successful supervisor, not just an output status file.
Original attempts/checkpoints are never changed. See docs/BRUNELLO_10K_SCORING.md.
"""
import argparse
from collections import Counter
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

PROJECT = Path('/moosefs/salehi/HPRCv2_CRISPRapido')
RUN = PROJECT/'runs/brunello_pam_aware_v1'
OLD = RUN/'01_pilot/timing-10000-v1/HG00735_pat'
PRIOR = PROJECT/'brunello_resume_20261123_v1'
OLD_PLAN = RUN/'00_manifest/timing-10000-execution-v3/execution-plan.json'
ROWS = 57578781
OLD_SIGNATURE = 'd5226c6c8002646581b135ebc4c1d13d9f5b9567dd26df6b5b5c4a57ca23fd9b'
SPOOL_SHA = 'b55abb94325d65af9eff565aea1d5ba2e626818e1115315e1579ef363f0b6615'
COLUMNS = 'row_index guide_id contig start_0based end_0based strand guide_sequence target_sequence pam cigar aligned_guide aligned_target mismatches bulge_events max_bulge_size cfd cfd_reason crispr_bulge crispr_bulge_reason'.split()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def publish(path, data):
    """Durable atomic no-clobber publication, including FAILED/CANDIDATE."""
    path = Path(path)
    tmp = path.with_name(path.name+'.tmp-'+uuid.uuid4().hex)
    with tmp.open('x') as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.link(tmp, path)
    tmp.unlink()
    fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def option(command, flag):
    require(command.count(flag) == 1, 'missing/repeated '+flag)
    return command[command.index(flag)+1]


def verify_files(files):
    for p, spec in files.items():
        require(Path(p).stat().st_size == spec['bytes'] and digest(p) == spec['sha256'], 'file drift: '+p)


def file_spec(path):
    return dict(bytes=Path(path).stat().st_size, sha256=digest(path))


def validate_plan(plan, root):
    require(plan['format'] == 'brunello-10k-scoring-migration-v1', 'wrong stage format')
    require(plan['assembly'] == 'HG00735_pat' and plan['guides'] == 10000 and plan['rows'] == ROWS, 'scope drift')
    require(plan['chunk_rows'] == 50000 and plan['concurrency'] == 1, 'unvalidated chunk/concurrency')
    require(plan['full_cohort_authorized'] is False and plan['search_allowed'] is False, 'search/cohort forbidden')
    require(plan['root'] == str(root) and plan['old_signature'] == OLD_SIGNATURE, 'stage identity drift')
    require(plan['old_plan'] == str(OLD_PLAN), 'wrong original plan')
    require(plan['input'] == str(OLD/'attempt-2874508/search/eligible.cfd.paf'), 'wrong source PAF')
    require(Path(plan['bundle']) == root/'scoring', 'wrong bundle')
    require(Path(plan['wrapper']) == root/'tooling/score_paf_chunked.py', 'unfrozen wrapper')
    require(Path(plan['spool']) == OLD/'attempt-2876332/score-stage/neural/scoring/records.jsonl.incomplete', 'wrong historical availability evidence')
    require(plan['spool_sha256'] == SPOOL_SHA, 'historical evidence pin drift')
    require(plan['original_outputs'][plan['input']]['sha256'] == 'fb09e1141f36e0f82561783dbacf6a302f92c2df547c2476a7ebe4edfb6da69c', 'original search digest drift')
    old = json.loads(OLD_PLAN.read_text())['commands']['score']; cmd = plan['command']
    require(cmd[1:4] == ['-I','-B',plan['wrapper']], 'not isolated frozen wrapper')
    for flag, value in {'--input':plan['input'], '--out-dir':plan['bundle'], '--scoring':'both', '--chunk-rows':'50000',
                        '--reference':option(old,'-r'), '--guides':option(old,'-G'),
                        '--model-root':option(old,'--crispr-bulge-model-root'), '--model-python':option(old,'--crispr-bulge-python'),
                        '--inner-script':option(old,'--scoring-script'), '--reference-index':option(old,'--scoring-reference-index'),
                        '--reference-index-sha256':option(old,'--scoring-reference-index-sha256'),
                        '--expected-reference-sha256':option(old,'--scoring-reference-sha256')}.items():
        require(option(cmd,flag) == value, 'scoring command drift: '+flag)
    require('--annotations' in cmd and '--preserve-cfd' in cmd, 'annotations/CFD must be retained')


def prepare(root, repository):
    require(root.parent == PROJECT and root.name.startswith('brunello_10k_scoring_'), 'unexpected stage root')
    root.mkdir(exist_ok=True)
    require(not (root/'plan.json').exists(), 'plan already frozen; do not overwrite')
    old = json.loads(OLD_PLAN.read_text())
    require(digest(OLD_PLAN) == OLD_SIGNATURE, 'old plan drift')
    require(old['guide_count'] == 10000 and old['assembly'] == 'HG00735_pat', 'not 10k pilot')
    prior_checks = json.loads((PRIOR/'verified-checkpoints.json').read_text())
    for name, spec in prior_checks.items():
        require(digest(OLD/'checkpoints'/(name+'.json')) == spec['sha256'], 'checkpoint metadata differs from reconciled pilot')
    source = json.loads((OLD/'checkpoints/search.json').read_text())
    audit = json.loads((OLD/'checkpoints/audit.json').read_text())
    require(source['signature'] == audit['signature'] == OLD_SIGNATURE and audit['rows'] == ROWS, 'checkpoint scope drift')
    tooling = root/'tooling'; tooling.mkdir()
    for src, name in [(repository/'research/crispr_bulge/score_paf_chunked.py', 'score_paf_chunked.py'),
                      (Path(__file__), 'scoring_completion.py'),
                      (Path(__file__).with_name('submit_scoring_completion.slurm'), 'submit.slurm'),
                      (PRIOR/'tooling/supervise.py', 'supervise.py')]:
        shutil.copyfile(src, tooling/name)
    require(digest(tooling/'score_paf_chunked.py') == digest(PRIOR/'tooling/score_paf_chunked.py'), 'wrapper differs from validated pilot')
    c = old['commands']['score']
    command = [sys.executable, '-I', '-B', str(tooling/'score_paf_chunked.py'),
               '--input', source['paf'], '--out-dir', str(root/'scoring'),
               '--scoring', 'both', '--annotations', '--preserve-cfd', '--chunk-rows', '50000',
               '--inner-script', option(c, '--scoring-script'), '--reference', option(c, '-r'),
               '--guides', option(c, '-G'), '--model-root', option(c, '--crispr-bulge-model-root'),
               '--model-python', option(c, '--crispr-bulge-python'),
               '--reference-index', option(c, '--scoring-reference-index'),
               '--reference-index-sha256', option(c, '--scoring-reference-index-sha256'),
               '--expected-reference-sha256', option(c, '--scoring-reference-sha256'),
               '--worker-timeout', '1800']
    pins = dict(old['pins'])
    pins[sys.executable] = digest(sys.executable)
    for p in tooling.iterdir():
        pins[str(p)] = digest(p)
    for name in ['pilot-validation.json', 'pilot-selection.json', 'injected-interruption.json']:
        pins[str(PRIOR/name)] = digest(PRIOR/name)
    plan = dict(format='brunello-10k-scoring-migration-v1', root=str(root), assembly='HG00735_pat',
                guides=10000, rows=ROWS, chunk_rows=50000, concurrency=1, full_cohort_authorized=False,
                search_allowed=False, authorization='User 2026-11-23: scoring-only completion of existing 10k pilot; preserve previous attempts; no SBDS/full cohort.',
                old_signature=OLD_SIGNATURE, old_plan=str(OLD_PLAN), input=source['paf'],
                bundle=str(root/'scoring'), wrapper=str(tooling/'score_paf_chunked.py'), command=command,
                checkpoints={str(OLD/'checkpoints'/name):file_spec(OLD/'checkpoints'/name) for name in ['search.json','audit.json']},
                original_outputs={**source['files'], **audit['files']}, pins=pins,
                spool=str(OLD/'attempt-2876332/score-stage/neural/scoring/records.jsonl.incomplete'),
                spool_sha256=SPOOL_SHA, spool_bytes=29881742568,
                rss_limit_kib=4194304, min_free_bytes=200_000_000_000,
                memory_request='8G', time_limit='36:00:00',
                note='New scoring identity only; old search/audit signatures untouched. Historical spool is conversion/availability evidence, never reused predictions.')
    validate_plan(plan, root)
    publish(root/'plan.json', plan)
    print(json.dumps(dict(plan=str(root/'plan.json'), sha256=digest(root/'plan.json'))))


def compare_bundles(a, b, tolerance=1e-6):
    """Exact rows/CFD/availability; only numeric CB may differ within tolerance."""
    peak = 0.; rows = 0
    for filename, cb_index in [('output.paf', None), ('annotations.tsv', 17)]:
        with (a/filename).open() as x, (b/filename).open() as y:
            count = 0
            for left in x:
                right = y.readline()
                require(bool(right), 'comparison row loss')
                lx, ry = left.rstrip('\n').split('\t'), right.rstrip('\n').split('\t')
                require(len(lx) == len(ry), 'comparison field loss')
                for i, (v, w) in enumerate(zip(lx, ry)):
                    numeric = (cb_index is None and v.startswith('cb:f:')) or (cb_index == i and v not in ('N/A', 'crispr_bulge'))
                    if numeric:
                        if cb_index is None:
                            require(w.startswith('cb:f:'), 'CB availability drift'); v, w = v[5:], w[5:]
                        u, z = float(v), float(w)
                        require(math.isfinite(u) and math.isfinite(z) and abs(u-z) <= tolerance, 'CB tolerance failure')
                        peak = max(peak, abs(u-z))
                    else:
                        require(v == w, 'non-CB field drift')
                count += 1
            require(not y.read(1), 'comparison extra rows')
            if filename == 'output.paf':
                rows = count
            else:
                require(count == rows+1, 'annotation/PAF comparison coverage')
    return dict(rows=rows, tolerance=tolerance, max_cb_absolute_delta=peak, status='PASS')


def check_prior_pilot():
    for name in ['single-pass', 'chunked']:
        d = PRIOR/name; m = json.loads((d/'metadata.json').read_text())
        require(digest(d/'output.paf') == m['output_sha256'], 'prior output hash drift')
        require(digest(d/'annotations.tsv') == m['annotations_sha256'], 'prior annotation hash drift')
    result = compare_bundles(PRIOR/'single-pass', PRIOR/'chunked')
    require(result['rows'] == 101000, 'prior coverage drift')
    with (PRIOR/'pilot-input.paf').open('rb') as f:
        for i, rows in enumerate([50000,50000,1000]):
            h = hashlib.sha256()
            for _ in range(rows):
                line = f.readline(); require(bool(line), 'prior input truncated'); h.update(line)
            bound = json.loads((PRIOR/f'chunked/chunk-config-{i:05d}.json').read_text())
            require(h.hexdigest() == bound['input_sha256'], 'prior chunk input hash drift')
        require(not f.read(1), 'prior input extra rows')
    result['checkpoint_input_hashes_verified'] = True
    return result


def verify_merged(source, output, annotations, spool, expected_rows, source_sha, spool_sha):
    """Independent full-stream check against source and historical conversion.

    The old failed run's complete conversion spool has no reusable predictions.
    Its archived digest and raw-row equality authenticate availability and every
    non-CB annotation field. Numeric CB predictions are newly computed.
    """
    counts = Counter(); hashes = [hashlib.sha256() for _ in range(4)]
    with open(source,'rb') as raw, open(output,'rb') as scored, open(annotations,'rb') as ann, open(spool,'rb') as historical:
        header = ann.readline(); hashes[2].update(header)
        require(header.decode().rstrip('\n').split('\t') == COLUMNS, 'annotation header drift')
        n = 0
        for line in raw:
            new, annotation, old = scored.readline(), ann.readline(), historical.readline()
            require(new and annotation and old, 'merged or historical row loss')
            for h, value in zip(hashes, [line,new,annotation,old]): h.update(value)
            record = json.loads(old)
            require(record['id'] == n and record['raw'].encode() == line, 'historical conversion order/input drift')
            prefix = line.rstrip(b'\r\n')+b'\t'
            require(new.startswith(prefix), 'PAF row/order/CFD changed')
            reason = record['status']; extra = new[len(prefix):].rstrip(b'\r\n').decode().split('\t')
            fields = annotation.decode().rstrip('\n').split('\t')
            require(len(fields) == len(COLUMNS), 'annotation field count')
            expected = [str(x) for x in record['annotation']]
            require(fields[:17] == expected[:17] and fields[18:] == expected[18:], 'non-CB annotation/CFD/availability drift')
            if reason == 'scored':
                require(len(extra) == 2 and extra[0].startswith('cb:f:') and extra[1] == 'bu:Z:scored', 'CB tag/availability drift')
                value = float(extra[0][5:])
                require(math.isfinite(value) and 0 <= value <= 1 and extra[0][5:] == fields[17], 'invalid CB score')
            else:
                require(extra == ['bu:Z:'+reason] and fields[17] == 'N/A', 'N/A must remain explicit')
            require(int(fields[0]) == n and fields[-1] == reason, 'annotation order/availability drift')
            counts[reason] += 1; n += 1
        require(not scored.read(1) and not ann.read(1) and not historical.read(1), 'extra merged/historical rows')
    require(n == expected_rows, 'final row coverage drift')
    require(hashes[0].hexdigest() == source_sha, 'source changed during scoring')
    require(hashes[3].hexdigest() == spool_sha, 'historical conversion digest drift')
    return dict(status='PASS', rows=n, availability=dict(counts), source_sha256=source_sha,
                output_sha256=hashes[1].hexdigest(), annotations_sha256=hashes[2].hexdigest(),
                historical_conversion_sha256=spool_sha, row_order_CFD_and_non_CB_annotations_unchanged=True)


def logical_bytes(directory):
    total = 0
    for p in directory.rglob('*'):
        try:
            if p.is_file(): total += p.stat().st_size
        except FileNotFoundError:
            pass  # wrapper-owned transient data may be cleaned concurrently
    return total


def work(root, plan, job):
    attempt = root/('attempt-'+job); attempt.mkdir()
    publish(attempt/'STARTED.json', dict(job=job, plan_sha256=digest(root/'plan.json'), start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    try:
        require(shutil.disk_usage(root).free >= plan['min_free_bytes'], 'output reserve preflight')
        verify_files(plan['checkpoints']); verify_files(plan['original_outputs'])
        for p, sha in plan['pins'].items(): require(digest(p) == sha, 'pin drift: '+p)
        require(Path(plan['spool']).stat().st_size == plan['spool_bytes'], 'historical evidence size drift')
        publish(attempt/'PREFLIGHT.json', dict(status='PASS', checkpoints_reused=True, prior_pilot=check_prior_pilot()))
        bundle = Path(plan['bundle'])
        # A previously completed wrapper may be revalidated after controller or
        # supervisor failure; it is not rescored and never trusted by status alone.
        complete = (bundle/'status').is_file() and (bundle/'status').read_text().startswith('COMPLETE')
        proc = None; start = time.monotonic(); last_storage = 0.; storage = logical_bytes(root); seen_chunks = set()
        with (attempt/'progress.jsonl').open('x', buffering=1) as log, (attempt/'wrapper.stdout').open('x') as out, (attempt/'wrapper.stderr').open('x') as err:
            if not complete:
                proc = subprocess.Popen(plan['command'], stdout=out, stderr=err)
            try:
                while proc is not None and proc.poll() is None:
                    free = shutil.disk_usage(root).free
                    require(free >= plan['min_free_bytes'], 'running output-reserve guard')
                    meta = {}; partial = bundle/'metadata.json.incomplete'
                    if partial.exists():
                        try: meta = json.loads(partial.read_text())
                        except (OSError, ValueError): pass
                    if time.monotonic()-last_storage >= 60:
                        storage = logical_bytes(root); last_storage = time.monotonic()
                    phase = 'scoring' if meta else 'duplicate_check' if (bundle/'chunks/anchor-check.sqlite.incomplete').exists() else 'preparing_or_merging'
                    new_metrics = [c for c in meta.get('chunks',[]) if c['index'] not in seen_chunks]
                    seen_chunks.update(c['index'] for c in new_metrics)
                    log.write(json.dumps(dict(new_chunk_metrics=new_metrics, utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()), elapsed_seconds=time.monotonic()-start,
                        phase=phase, chunks_complete=len(meta.get('chunks',[])), rows_complete=sum(c['rows'] for c in meta.get('chunks',[])),
                        expected_chunks=1152, logical_stage_bytes=storage, free_bytes=free))+'\n')
                    time.sleep(10)
                require(proc is None or proc.returncode == 0, 'bounded wrapper failed; retained outputs/logs')
            finally:
                if proc is not None and proc.poll() is None:
                    proc.terminate()
                    try: proc.wait(timeout=5)
                    except subprocess.TimeoutExpired: proc.kill(); proc.wait()
        meta = json.loads((bundle/'metadata.json').read_text())
        require((bundle/'status').read_text().startswith('COMPLETE') and meta['rows'] == ROWS, 'wrapper incomplete/coverage drift')
        publish(attempt/'SCORING_FINISHED.json', dict(chunks=len(meta['chunks']), rows=meta['rows'], wrapper_wall_this_attempt=time.monotonic()-start,
            logical_stage_bytes=logical_bytes(root), storage_caveat='logical bytes, not physical replication or allocation'))
        source_sha = plan['original_outputs'][plan['input']]['sha256']
        result = verify_merged(plan['input'], bundle/'output.paf', bundle/'annotations.tsv', plan['spool'], ROWS, source_sha, plan['spool_sha256'])
        require(result['output_sha256'] == meta['output_sha256'] and result['annotations_sha256'] == meta['annotations_sha256'], 'merged hashes disagree')
        verify_files(plan['checkpoints'])
        for p, sha in plan['pins'].items(): require(digest(p) == sha, 'post-score pin drift: '+p)
        publish(attempt/'VALIDATION.json', result)
        files = {str(bundle/name):file_spec(bundle/name) for name in ['metadata.json','status']}
        files[str(bundle/'output.paf')] = dict(bytes=(bundle/'output.paf').stat().st_size, sha256=result['output_sha256'])
        files[str(bundle/'annotations.tsv')] = dict(bytes=(bundle/'annotations.tsv').stat().st_size, sha256=result['annotations_sha256'])
        files[str(attempt/'VALIDATION.json')] = file_spec(attempt/'VALIDATION.json')
        publish(attempt/'CANDIDATE.json',dict(status='AWAITING_SUPERVISOR_AND_AFTEROK', job=job,
            plan_sha256=digest(root/'plan.json'), validation=result, files=files))
    except BaseException as exc:
        publish(attempt/'FAILED.json', dict(error=repr(exc), job=job))
        raise


def finalize(root, plan, score_job):
    require(os.environ.get('SLURM_JOB_ID') != score_job, 'finalize requires separate afterok job')
    state = subprocess.check_output(['scontrol','show','job',score_job], text=True)
    require('JobState=COMPLETED' in state and 'ExitCode=0:0' in state, 'scoring Slurm job not completed successfully')
    report = root/('supervisor-'+score_job+'.json'); supervisor = json.loads(report.read_text())
    require(supervisor['reason'] == 'completed' and supervisor['child_returncode'] == 0, 'supervisor failed')
    candidate = json.loads((root/('attempt-'+score_job)/'CANDIDATE.json').read_text())
    require(candidate['job'] == score_job and candidate['plan_sha256'] == digest(root/'plan.json'), 'candidate identity drift')
    verify_files(candidate['files']); verify_files(plan['checkpoints'])
    publish(root/'SUCCESS.json',dict(status='COMPLETE', scoring_job=score_job, finalizer_job=os.environ['SLURM_JOB_ID'],
        plan_sha256=digest(root/'plan.json'), supervisor=supervisor, supervisor_sha256=digest(report),
        scheduler_evidence=state, validation=candidate['validation'], files=candidate['files'],
        scope='Existing 10k HG00735_pat scoring only; no SBDS/search/cohort run'))
    (root/'RESULTS.md').write_text('# Brunello 10k scoring COMPLETE\n\n'+
        f'Scoring job {score_job} COMPLETED 0:0; afterok finalizer {os.environ["SLURM_JOB_ID"]}.\n\n'+
        f'Verified {ROWS:,} rows: source order/coverage/CFD and historical CB availability unchanged.\n\n'+
        'See SUCCESS.json, attempt VALIDATION.json, chunk metadata and supervisor/progress logs. Original attempts/checkpoints preserved. No SBDS or full cohort.\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['prepare','work','finalize'])
    p.add_argument('--root',type=Path,required=True); p.add_argument('--repository',type=Path)
    p.add_argument('--score-job'); a = p.parse_args()
    if a.mode == 'prepare':
        require(a.repository is not None, 'repository required'); prepare(a.root,a.repository); return
    plan = json.loads((a.root/'plan.json').read_text()); validate_plan(plan,a.root)
    require(os.environ.get('SLURM_CPUS_PER_TASK') == '1' and os.environ.get('SLURM_JOB_ID'), 'one-CPU Slurm job required')
    with (a.root/'stage.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if a.mode == 'work': work(a.root,plan,os.environ['SLURM_JOB_ID'])
        else: require(bool(a.score_job),'score job required'); finalize(a.root,plan,a.score_job)


if __name__ == '__main__':
    main()
