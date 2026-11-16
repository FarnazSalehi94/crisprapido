#!/usr/bin/env python3
"""Bounded 10k Brunello worker; never submits and cannot run a full cohort.

Immutable attempts + independent search/audit/score checkpoints. A successful
worker writes a candidate only; finalize requires successful external supervisor
evidence. All published checkpoints are exclusive, durable, and hash-verified.
"""
import argparse
import collections
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import tempfile
import time
import uuid

ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/brunello_pam_aware_v1')
CONFIG = dict(m=4, b=1, z=1, f=.75, p='GG')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def publish(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    with temporary.open('x') as f:
        json.dump(value, f, indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.link(temporary, path)  # atomic no-clobber publication
    temporary.unlink()  # only this invocation's metadata temporary
    fd = os.open(path.parent, os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)


def option(command, flag):
    require(command.count(flag) == 1, 'missing/repeated flag ' + flag)
    return command[command.index(flag)+1]


def replace(command, **flags):
    command = list(command)
    for flag, value in flags.items():
        option(command, flag)
        command[command.index(flag)+1] = str(value)
    return command


def verify_files(files):
    require(bool(files), 'empty file manifest')
    for path, spec in files.items():
        require(Path(path).is_file() and not Path(path).is_symlink(), 'nonregular output: '+path)
        require(Path(path).stat().st_size == spec['bytes'] and digest(path) == spec['sha256'], 'output drift: '+path)


def inventory(directory):
    files = {}
    for path in sorted(directory.rglob('*')):
        if path.is_file():
            require(not path.is_symlink(), 'symlink in output inventory')
            with path.open('rb') as f: os.fsync(f.fileno())
            files[str(path)] = dict(bytes=path.stat().st_size, sha256=digest(path))
    require(bool(files), 'empty stage output')
    return files


def checkpoint(path, signature, execute):
    if path.exists():
        data = json.loads(path.read_text())
        require(data['signature'] == signature, 'checkpoint signature drift')
        verify_files(data['files'])
        return data
    data = execute()
    data['signature'] = signature
    verify_files(data['files'])
    publish(path, data)
    return data


def verify_pins(plan):
    for path, sha in plan['pins'].items():
        require(digest(path) == sha, 'input/tool pin drift: '+path)


def validate_plan(plan):
    require(plan['format'] == 'brunello-timing-10000-v1', 'not bounded timing plan')
    require(plan['guide_count'] == 10000 and plan['assembly'] == 'HG00735_pat', 'benchmark scope drift')
    require(plan['configuration'] == CONFIG and plan['full_cohort_authorized'] is False, 'configuration/cohort drift')
    require(bool(plan['authorization']), 'benchmark authorization required')
    require(Path(plan['output']) == ROOT/'01_pilot/timing-10000-v1/HG00735_pat', 'cross-project output')
    c = plan['commands']; s = c['search']
    for flag, value in {'-m':'4','-b':'1','-z':'1','-f':'0.75','-p':'GG','-t':'1','--search-backend':'pam-anchored'}.items():
        require(option(s, flag) == value, 'search configuration drift: '+flag)
    ref, guides = option(s, '-r'), option(s, '-G')
    require(option(c['score'], '--scoring') == 'both', 'both scores required')
    for command, rf, gf in [(c['score'], '-r', '-G'), (c['audit'], '--fasta', '--guides')]:
        require(option(command, rf) == ref and option(command, gf) == guides, 'stage input drift')
    for flag, value in {'--max-mismatches':'4','--max-bulges':'1','--max-bulge-size':'1'}.items():
        require(option(c['audit'], flag) == value, 'audit domain drift')
    for path in (s[0], option(s, '--pam-search-bin'), ref, ref+'.fai', guides,
                 c['audit'][1], option(c['audit'], '--samtools'), option(c['score'], '--scoring-script')):
        require(path in plan['pins'], 'missing required pin: '+path)
    lines = Path(guides).read_text().splitlines()
    require(len(lines) == 20000, 'guide count')
    names, seqs = lines[::2], lines[1::2]
    require(len(set(names)) == len(set(seqs)) == 10000, 'duplicate guide identity')
    require(all(n.startswith('>') for n in names) and all(len(s)==20 and not set(s)-set('ACGT') for s in seqs), 'guide alphabet/format')


def command_run(command, directory, filename):
    directory.mkdir(exist_ok=False)
    before = resource.getrusage(resource.RUSAGE_CHILDREN); start = time.monotonic()
    output = directory/filename
    with output.open('x') as out, (directory/'stderr.log').open('x') as err:
        subprocess.run(list(map(str, command)), stdout=out, stderr=err, check=True, timeout=60*60*60)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    publish(directory/'metrics.json', dict(command=list(map(str, command)), wall_seconds=time.monotonic()-start,
        user_seconds=after.ru_utime-before.ru_utime, system_seconds=after.ru_stime-before.ru_stime,
        cumulative_max_child_rss_kib=after.ru_maxrss, rss_caveat='cumulative individual descendant maximum; tree RSS comes from supervisor'))
    return output


def score_counts(annotations, expected):
    counts = collections.Counter()
    with annotations.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            require(row['pam'][1:] == 'GG', 'score PAM drift')
            require(row['cfd'] != 'N/A' and math.isfinite(float(row['cfd'])) and 0 <= float(row['cfd']) <= 1, 'invalid CFD')
            reason = row['crispr_bulge_reason']
            require(bool(reason), 'missing score availability reason')
            if reason == 'scored':
                require(row['crispr_bulge'] != 'N/A' and math.isfinite(float(row['crispr_bulge'])) and 0 <= float(row['crispr_bulge']) <= 1, 'invalid defined CB score')
            else:
                require(row['crispr_bulge'] == 'N/A', 'unsupported score must be explicit N/A')
            counts['rows'] += 1
            counts['cb_'+reason] += 1
    require(counts['rows'] == expected, 'scoring row loss')
    return dict(counts)


def verify_fai(plan, scratch):
    ref = Path(option(plan['commands']['search'], '-r'))
    alias = scratch/'reference.fa'; alias.symlink_to(ref)
    subprocess.run([option(plan['commands']['audit'], '--samtools'), 'faidx', str(alias)], check=True, timeout=900)
    require(digest(str(alias)+'.fai') == plan['pins'][str(ref)+'.fai'], 'FAI/reference provenance mismatch')
    return alias


def work(plan, signature, job, scratch_base=None):
    out = Path(plan['output']); stages = out/'checkpoints'; stages.mkdir(exist_ok=True)
    success = out/'SUCCESS.json'
    if success.exists():
        data = json.loads(success.read_text()); require(data['signature'] == signature, 'SUCCESS drift')
        verify_files(data['files']); return
    attempt = out/('attempt-'+job); attempt.mkdir(exist_ok=False)
    publish(attempt/'STARTED.json', dict(signature=signature, job=job))
    scratch = None
    try:
        base = Path(scratch_base) if scratch_base else Path('/scratch')/os.environ['USER']
        require(shutil.disk_usage(base.parent).free >= 60_000_000_000, 'scratch preflight')
        require(shutil.disk_usage(out).free >= 150_000_000_000, 'output reserve preflight')
        base.mkdir(exist_ok=True)
        scratch = Path(tempfile.mkdtemp(prefix='brunello10k-'+job+'-', dir=base))
        os.environ.update(TMPDIR=str(scratch), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                          TF_NUM_INTRAOP_THREADS='1', TF_NUM_INTEROP_THREADS='1', CUDA_VISIBLE_DEVICES='-1')
        c = plan['commands']; ref = Path(option(c['search'], '-r'))
        # Recompute .fai from the pinned reference; never write beside original.
        alias = verify_fai(plan, scratch)
        publish(attempt/'fai-validation.json', dict(reference=str(ref), fai_sha256=digest(str(alias)+'.fai')))

        def search():
            paf = command_run(c['search'], attempt/'search', 'eligible.cfd.paf')
            return dict(paf=str(paf), files=inventory(paf.parent))
        search_data = checkpoint(stages/'search.json', signature, search)
        paf = Path(search_data['paf'])

        def audit():
            directory = attempt/'audit-stage'
            cmd = replace(c['audit'], **{'--paf':paf, '--fasta':alias, '--output':directory/'audit'})
            command_run(cmd, directory, 'stdout.json')
            summary = json.loads((directory/'audit/summary.json').read_text())
            require(summary['status'] == 'PASS' and summary['errors'] == 0 and (directory/'audit/PASS').is_file(), 'audit failed')
            require(summary['counts']['verified'] == summary['rows'], 'audit row accounting')
            return dict(rows=summary['rows'], files=inventory(directory))
        audit_data = checkpoint(stages/'audit.json', signature, audit)

        def score():
            directory = attempt/'score-stage'
            cmd = replace(c['score'], **{'--score-paf':paf, '--scoring-output-dir':directory/'neural'})
            # Execution-layer override only: same plan identity, so search/audit
            # checkpoints stay valid. The replacement scorer must wrap the same
            # pinned inner implementation (chunked wrapper); the actual command
            # is recorded in metrics.json for provenance.
            override = os.environ.get('BRUNELLO_SCORING_SCRIPT')
            override_sha = os.environ.get('BRUNELLO_SCORING_SCRIPT_SHA256')
            if override:
                require(Path(override).is_file() and not Path(override).is_symlink(),
                        'scoring override missing: '+override)
                # The override must be the reviewed chunked wrapper exactly;
                # a digest mismatch refuses to run rather than changing science.
                require(override_sha and digest(override) == override_sha,
                        'scoring override digest mismatch: '+override)
            scored = command_run(cmd, directory, 'eligible.both.paf')
            require(digest(scored) == digest(directory/'neural/scoring/output.paf'), 'scoring stdout/bundle mismatch')
            counts = score_counts(directory/'neural/scoring/annotations.tsv', audit_data['rows'])
            with scored.open() as f:
                require(sum(1 for _ in f) == audit_data['rows'], 'scored PAF row loss')
            metadata = json.loads((directory/'neural/scoring/metadata.json').read_text())
            require(metadata['input_sha256'] == digest(paf), 'score input identity')
            return dict(counts=counts, files=inventory(directory))
        score_data = checkpoint(stages/'score.json', signature, score)
        verify_pins(plan)
        files = {}
        for data in (search_data, audit_data, score_data): files.update(data['files'])
        files.update(inventory(stages))
        files.update({str(attempt/'fai-validation.json'):dict(bytes=(attempt/'fai-validation.json').stat().st_size,
                    sha256=digest(attempt/'fai-validation.json'))})
        publish(attempt/'CANDIDATE.json', dict(status='AWAITING_SUPERVISOR', signature=signature,
            job=job, files=files, counts=score_data['counts']))
        shutil.rmtree(scratch)  # only successful job-private scratch
    except BaseException as exc:
        publish(attempt/'FAILED.json', dict(error=str(exc), scratch=str(scratch), signature=signature))
        raise


def finalize(out, signature, job, supervisor):
    success = out/'SUCCESS.json'
    if success.exists():
        data = json.loads(success.read_text()); require(data['signature'] == signature, 'SUCCESS drift')
        verify_files(data['files']); return
    data = json.loads((out/('attempt-'+job)/'CANDIDATE.json').read_text())
    require(data['signature'] == signature and data['job'] == job, 'candidate identity drift')
    report = json.loads(supervisor.read_text())
    require(report['reason'] == 'completed' and report['child_returncode'] == 0, 'supervisor did not succeed')
    require(supervisor == out/('supervisor-'+job+'.json'), 'supervisor path/job drift')
    verify_files(data['files'])
    data.update(status='COMPLETE', supervisor=str(supervisor), supervisor_sha256=digest(supervisor))
    data['files'][str(supervisor)] = dict(bytes=supervisor.stat().st_size, sha256=digest(supervisor))
    publish(success, data)


def main():
    p = argparse.ArgumentParser(); p.add_argument('--plan', type=Path, required=True)
    modes = p.add_mutually_exclusive_group()
    modes.add_argument('--finalize', type=Path)
    modes.add_argument('--preflight-only', action='store_true'); a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_CPUS_PER_TASK') == '1', 'one-CPU Slurm required')
    plan = json.loads(a.plan.read_text()); signature = digest(a.plan)
    validate_plan(plan)
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=True)
    with (out/'benchmark.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verify_pins(plan)
        if a.preflight_only:
            require(shutil.disk_usage('/scratch').free >= 60_000_000_000, 'scratch preflight')
            require(shutil.disk_usage(out).free >= 150_000_000_000, 'output reserve preflight')
            base = Path('/scratch')/os.environ['USER']; base.mkdir(exist_ok=True)
            with tempfile.TemporaryDirectory(prefix='brunello10k-preflight-', dir=base) as tmp:
                alias = verify_fai(plan, Path(tmp))
                interpreter = Path(plan.get('recorded_not_pinned') or plan['commands']['audit'][0])
                publish(out/('preflight-'+os.environ['SLURM_JOB_ID']+'.json'),
                    dict(status='PASS', signature=signature, all_input_pins_verified=True,
                         audit_interpreter=str(interpreter), audit_interpreter_sha256=digest(interpreter),
                         fai_sha256=digest(str(alias)+'.fai'), scratch_free=shutil.disk_usage('/scratch').free,
                         output_free=shutil.disk_usage(out).free))
        elif a.finalize: finalize(out, signature, os.environ['SLURM_JOB_ID'], a.finalize)
        else: work(plan, signature, os.environ['SLURM_JOB_ID'])


if __name__ == '__main__':
    main()
