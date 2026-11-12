import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

WORKER = Path(__file__).with_name('cohort_worker.py')

FAKE_CLI = r'''#!/usr/bin/env python3
import os, sys
args = sys.argv[1:]
CLASS = {'GG': 'NGG', 'AG': 'NAG', 'GA': 'NGA', 'AA': 'NAA'}
if '--score-paf' in args:
    out = args[args.index('--scoring-output-dir')+1]
    os.makedirs(out+'/scoring', exist_ok=True)
    src = args[args.index('--score-paf')+1]
    pam_class = next(p for p in ('NGG','NAG','NGA','NAA') if '/'+p+'/' in src)
    suffix = dict(zip(('NGG','NAG','NGA','NAA'), ('GG','AG','GA','AA')))[pam_class]
    rows = open(src).read().splitlines()
    open(out+'/scoring/output.paf','w').write(''.join(r+'\n' for r in rows))
    with open(out+'/scoring/annotations.tsv','w') as f:
        f.write('row_index\tguide_id\tcontig\tstrand\tstart_0based\tend_0based\tpam\tcfd\tcrispr_bulge\tcrispr_bulge_reason\n')
        for i, r in enumerate(rows):
            f.write(f'{i}\tg\tc\t+\t0\t{20+i}\tA{suffix}\t0.5\tN/A\talignment_exceeds_24_columns\n')
    sys.stdout.write(''.join(r+'\n' for r in rows))
else:
    suffix = args[args.index('-p')+1]
    pam_class = CLASS[suffix]
    if os.environ.get('FAKE_FAIL_PAM') == pam_class:
        sys.exit(3)
    for i in range(3):
        print(f'g\t20\t0\t20\t+\tc\t100\t{5+i}\t{25+i}\t20\t20\t255\tas:i:0\tnm:i:0\tng:i:0\tbs:i:0\tcg:Z:20=\tcf:f:1.0000')
'''

FAKE_AUDITOR = r'''
def audit(paf, fasta, guides, matrix_dir, samtools, output, **kw):
    from pathlib import Path
    import json
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    rows = sum(1 for _ in open(paf))
    (output/'records.tsv').write_text('x\n'*rows)
    summary = dict(rows=rows, status='PASS')
    (output/'summary.json').write_text(json.dumps(summary))
    return summary
'''


def digest(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(65536), b''): h.update(b)
    return h.hexdigest()


class CohortWorkerTests(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.root = Path(self.t.name)/'sbds_ngg_nrr_v1'
        for d in ('02_search', '03_scoring', '04_qc'):
            (self.root/d).mkdir(parents=True)
        self.cli = self.root/'fake-cli'; self.cli.write_text(FAKE_CLI); self.cli.chmod(0o700)
        self.auditor = self.root/'fake_auditor.py'; self.auditor.write_text(FAKE_AUDITOR)
        self.guides = self.root/'guides.fa'; self.guides.write_text('>SBDS\nACTGAAATCTGTAAGCAGGC\n')
        self.source = self.root/'source'; self.source.mkdir()
        for n in ('mismatch_scores.txt', 'pam_scores.txt'):
            (self.source/n).write_text('x')
        self.ref = self.root/'reference.fa'; self.ref.write_text('>c\n'+'ACGT'*10+'\n')
        self.fai = self.root/'reference.fa.fai'; self.fai.write_text('c\t40\n')
        (self.root/'index.crispr-binding.json').write_text('{}')
        ri = self.root/'ri.json'; ri.write_text('{}')
        binding = self.root/'index.crispr-binding.json'
        self.plan = dict(format='sbds-nrr-cohort-v1', root=str(self.root),
            source=str(self.source), guides=str(self.guides), models='/m', model_python=sys.executable,
            samtools='/usr/bin/true', auditor=str(self.auditor),
            assemblies=dict(TEST=dict(assembly='TEST', cli=str(self.cli), native=str(self.cli),
                index=str(self.root/'index'), reference=str(self.ref), reference_sha256=digest(self.ref),
                fai=str(self.fai), reference_index=str(ri),
                pins={str(self.ref): digest(self.ref), str(binding): digest(binding), str(ri): digest(ri)})),
            pins={str(p): digest(p) for p in [self.cli, self.auditor, self.guides,
                self.source/'mismatch_scores.txt', self.source/'pam_scores.txt']},
            limits=dict(rss_kib=1048576, min_start_available_kib=1, min_running_available_kib=1,
                min_scratch_bytes=1, min_running_disk_bytes=1, stage_timeout_seconds=60, file_limit_bytes=10**9))
        self.plan_file = self.root/'plan.json'
        self.plan_file.write_text(json.dumps(self.plan))
        self.scratch = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.scratch, ignore_errors=True)
        self.t.cleanup()

    def env(self, job='1'):
        e = dict(os.environ, SLURM_JOB_ID=job, SLURM_CPUS_PER_TASK='1', SLURMD_NODENAME='x')
        for k in list(e):
            if k.startswith('FAKE_'): del e[k]
        return e

    def invoke(self, job, env=None):
        return subprocess.run([sys.executable, str(WORKER), '--plan', str(self.plan_file), '--assembly', 'TEST'],
                               capture_output=True, text=True, env=dict(env or self.env(job), TMPDIR=str(self.scratch)))

    def test_first_pam_failure_then_clean_restart(self):
        r = self.invoke('1', dict(self.env('1'), FAKE_FAIL_PAM='NGG'))
        self.assertNotEqual(r.returncode, 0)
        failed = json.loads((self.root/'02_search/TEST/attempt-1/FAILED.json').read_text())
        self.assertIn('stage returned 3', failed['error'])
        self.assertFalse((self.root/'02_search/TEST/NGG/SEARCH_COMPLETE.json').exists())
        r2 = self.invoke('2')
        self.assertEqual(r2.returncode, 0, r2.stderr)
        s = json.loads((self.root/'02_search/TEST/SUCCESS.json').read_text())
        self.assertEqual(s['status'], 'COMPLETE')
        self.assertEqual(len(s['outputs']), len({k: v['sha256'] for k, v in s['outputs'].items()}))
        r3 = self.invoke('3')
        self.assertEqual(r3.returncode, 0, r3.stderr)
        self.assertIn('Already complete', r3.stdout)

    def test_midrun_failure_keeps_completed_pam_checkpoints(self):
        # unmodified FAKE_CLI with FAKE_FAIL_PAM=NGA fails only the third PAM's
        # search; NGG/NAG must complete and be reused on restart.
        r = self.invoke('1', dict(self.env('1'), FAKE_FAIL_PAM='NGA'))
        self.assertNotEqual(r.returncode, 0)
        self.assertTrue((self.root/'02_search/TEST/NGG/SEARCH_COMPLETE.json').exists())
        self.assertTrue((self.root/'02_search/TEST/NAG/SEARCH_COMPLETE.json').exists())
        self.assertFalse((self.root/'02_search/TEST/NGA/SEARCH_COMPLETE.json').exists())
        # uncommitted partial NGA output exists without checkpoint -> removed on rerun
        (self.root/'02_search/TEST/NGA/eligible.cfd.paf').write_text('garbage\n')
        self.cli.write_text(FAKE_CLI)
        self.plan['pins'][str(self.cli)] = digest(self.cli)
        self.plan_file.write_text(json.dumps(self.plan))
        r2 = self.invoke('2')
        self.assertEqual(r2.returncode, 0, r2.stderr)
        reused = [json.loads((self.root/'02_search/TEST/attempt-2'/f'{p}.search-metrics.json').read_text())['reused_search']
                  for p in ('NGG', 'NAG')]
        self.assertEqual(reused, [True, True])
        s = json.loads((self.root/'02_search/TEST/SUCCESS.json').read_text())
        self.assertTrue(all(k in s['outputs'] for k in
            [str(self.root/'02_search/TEST'/p/'eligible.cfd.paf') for p in ('NGG','NAG','NGA','NAA')]))

    def test_tampered_success_output_is_detected(self):
        self.assertEqual(self.invoke('1').returncode, 0)
        target = next(p for p in json.loads((self.root/'02_search/TEST/SUCCESS.json').read_text())['outputs']
                     if p.endswith('eligible.cfd.paf'))
        Path(target).write_text('tampered\n')
        r = self.invoke('2')
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('completed output changed', r.stderr)

    def test_stage_timeout_kills_process_group(self):
        slow = self.root/'slow-cli'
        slow.write_text('#!/usr/bin/env python3\nimport time; time.sleep(120)\n'); slow.chmod(0o700)
        self.plan['assemblies']['TEST']['cli'] = str(slow)
        self.plan['pins'][str(slow)] = digest(slow)
        self.plan['limits']['stage_timeout_seconds'] = 2
        self.plan_file.write_text(json.dumps(self.plan))
        r = self.invoke('1')
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('stage timeout', json.loads((self.root/'02_search/TEST/attempt-1/FAILED.json').read_text())['error'])
        out = subprocess.run(['pgrep', '-f', str(slow)], capture_output=True, text=True)
        self.assertEqual(out.stdout, '', 'slow process leaked')

    def test_bad_plan_or_binding_rejected(self):
        bad = dict(self.plan); bad['format'] = 'other'
        pf = self.root/'bad.json'; pf.write_text(json.dumps(bad))
        r = subprocess.run([sys.executable, str(WORKER), '--plan', str(pf), '--assembly', 'TEST'],
                           capture_output=True, text=True, env=self.env('1'))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('unsupported plan', r.stderr)
        os.rename(self.root/'index.crispr-binding.json', self.root/'binding-moved.json')
        r2 = self.invoke('1')
        self.assertNotEqual(r2.returncode, 0)
        self.assertIn('binding sidecar missing', r2.stderr)

    def test_second_job_same_assembly_waits_or_fails_not_duplicates(self):
        import fcntl
        (self.root/'02_search/TEST').mkdir(parents=True, exist_ok=True)
        lock = (self.root/'02_search/TEST/assembly.lock').open('a')
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        r = self.invoke('1')
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('resource temporarily unavailable', r.stderr.lower())
        fcntl.flock(lock, fcntl.LOCK_UN)


if __name__ == '__main__':
    unittest.main()
