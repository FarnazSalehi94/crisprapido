import fcntl
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('timing', Path(__file__).with_name('timing_worker.py'))
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)


class CheckpointTests(unittest.TestCase):
    def test_unvalidated_scoring_override_rejected(self):
        with patch.dict(os.environ, {'BRUNELLO_SCORING_SCRIPT': '/some/wrapper'}):
            with self.assertRaisesRegex(ValueError, 'override disabled'):
                w.reject_scoring_override()
        with patch.dict(os.environ, {}, clear=True):
            w.reject_scoring_override()

    def test_reuse_tampering_and_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); output = root/'stage'; output.mkdir(); (output/'paf').write_text('hit\n')
            calls = []
            def execute():
                calls.append(1); return dict(files=w.inventory(output))
            cp = root/'checkpoint.json'
            w.checkpoint(cp, 'sig', execute); w.checkpoint(cp, 'sig', execute)
            self.assertEqual(len(calls), 1)
            with self.assertRaises(ValueError): w.checkpoint(cp, 'changed', execute)
            (output/'paf').write_text('tampered')
            with self.assertRaises(ValueError): w.checkpoint(cp, 'sig', execute)

    def test_crash_before_publication_retries(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); stage = root/'stage'; stage.mkdir(); (stage/'paf').write_text('hit\n')
            cp = root/'cp.json'
            with patch.object(w.os, 'link', side_effect=OSError('injected crash')):
                with self.assertRaises(OSError): w.checkpoint(cp, 'sig', lambda: dict(files=w.inventory(stage)))
            self.assertFalse(cp.exists())
            w.checkpoint(cp, 'sig', lambda: dict(files=w.inventory(stage)))
            self.assertTrue(cp.exists())

    def test_failed_command_no_checkpoint(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            def execute():
                w.command_run([sys.executable, '-c', 'print("partial"); raise SystemExit(9)'], root/'stage', 'paf')
            with self.assertRaises(subprocess.CalledProcessError): w.checkpoint(root/'cp.json', 'sig', execute)
            self.assertFalse((root/'cp.json').exists()); self.assertTrue((root/'stage/paf').exists())

    def test_pin_drift(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'tool'; p.write_text('good'); plan = dict(pins={str(p):w.digest(p)})
            w.verify_pins(plan); p.write_text('changed')
            with self.assertRaises(ValueError): w.verify_pins(plan)

    def test_publication_never_overwrites(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'cp'; w.publish(p, {'x':1})
            with self.assertRaises(FileExistsError): w.publish(p, {'x':2})
            self.assertEqual(json.loads(p.read_text()), {'x':1})

    def test_score_availability_and_loss(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'annotations.tsv'
            p.write_text('pam\tcfd\tcrispr_bulge\tcrispr_bulge_reason\nAGG\t0.5\t0.2\tscored\nTGG\t0.1\tN/A\tunsupported\n')
            self.assertEqual(w.score_counts(p, 2)['cb_unsupported'], 1)
            with self.assertRaises(ValueError): w.score_counts(p, 3)
            p.write_text(p.read_text().replace('0.2\tscored', 'N/A\tscored'))
            with self.assertRaises(ValueError): w.score_counts(p, 2)


class WorkflowTests(unittest.TestCase):
    def fixture(self, root):
        out = root/'out'; out.mkdir(); ref = root/'reference.fa'; ref.write_text('>chr\nACGT\n')
        fai = Path(str(ref)+'.fai'); fai.write_text('chr\t4\t5\t4\t5\n')
        return dict(output=str(out), pins={str(fai):w.digest(fai)}, commands={
            'search':['search','-r',str(ref)],
            'audit':['audit','--samtools','samtools','--paf','old','--fasta',str(ref),'--output','old'],
            'score':['score','--score-paf','old','--scoring-output-dir','old']})

    def simulate(self, plan, calls, fail_stage=None):
        def run(command, directory, filename):
            stage = command[0]; calls.append(stage); directory.mkdir()
            out = directory/filename; out.write_text('hit\n')
            if stage == fail_stage: raise RuntimeError('injected '+stage+' failure')
            if stage == 'audit':
                dest = Path(w.option(command, '--output')); dest.mkdir()
                (dest/'summary.json').write_text(json.dumps(dict(status='PASS', errors=0, rows=1, counts=dict(verified=1))))
                (dest/'PASS').write_text('PASS')
            if stage == 'score':
                dest = Path(w.option(command, '--scoring-output-dir'))/'scoring'; dest.mkdir(parents=True)
                (dest/'output.paf').write_text('hit\n')
                (dest/'annotations.tsv').write_text('pam\tcfd\tcrispr_bulge\tcrispr_bulge_reason\nAGG\t0.5\t0.2\tscored\n')
                (dest/'metadata.json').write_text(json.dumps(dict(input_sha256=w.digest(w.option(command, '--score-paf')))))
            return out
        return run

    def fake_faidx(self, command, **kwargs):
        alias = Path(command[-1]); Path(str(alias)+'.fai').write_text('chr\t4\t5\t4\t5\n')

    def run_work(self, plan, root, calls, job, fail=None):
        with patch.object(w, 'command_run', side_effect=self.simulate(plan, calls, fail)), patch.object(w.subprocess, 'run', side_effect=self.fake_faidx):
            with patch.object(w.shutil, 'disk_usage') as disk:
                disk.return_value.free = 10**13
                w.work(plan, 'sig', job, scratch_base=root/'scratch')

    def test_score_failure_resume_and_supervisor_gate(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); plan = self.fixture(root); calls = []
            with self.assertRaises(RuntimeError): self.run_work(plan, root, calls, '1', 'score')
            self.assertFalse((root/'out/SUCCESS.json').exists())
            self.run_work(plan, root, calls, '2')
            self.assertEqual(calls, ['search','audit','score','score'])
            out = root/'out'; supervisor = out/'supervisor-2.json'
            supervisor.write_text(json.dumps(dict(reason='sampled_rss_limit', child_returncode=0)))
            with self.assertRaises(ValueError): w.finalize(out, 'sig', '2', supervisor)
            self.assertFalse((out/'SUCCESS.json').exists())
            supervisor.write_text(json.dumps(dict(reason='completed', child_returncode=0)))
            w.finalize(out, 'sig', '2', supervisor)
            self.assertEqual(json.loads((out/'SUCCESS.json').read_text())['status'], 'COMPLETE')
            w.work(plan, 'sig', '3', scratch_base=root/'scratch')
            self.assertFalse((out/'attempt-3').exists())
            # Completed outputs remain checksummed on skip.
            (out/'attempt-2/score-stage/eligible.both.paf').write_text('changed')
            with self.assertRaises(ValueError): w.work(plan, 'sig', '4', scratch_base=root/'scratch')

    def test_audit_failure_reuses_search(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); plan = self.fixture(root); calls = []
            with self.assertRaises(RuntimeError): self.run_work(plan, root, calls, '1', 'audit')
            self.run_work(plan, root, calls, '2')
            self.assertEqual(calls, ['search','audit','audit','score'])

    def test_fai_drift_blocks_search(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); plan = self.fixture(root); calls = []
            plan['pins'][next(iter(plan['pins']))] = '0'*64
            with self.assertRaises(ValueError): self.run_work(plan, root, calls, '1')
            self.assertEqual(calls, [])

    def test_lock_contention(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); plan = self.fixture(root); p = root/'plan.json'; p.write_text(json.dumps(plan))
            with (root/'out/benchmark.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                with patch.dict(os.environ, SLURM_JOB_ID='test', SLURM_CPUS_PER_TASK='1'), patch.object(sys, 'argv', ['worker','--plan',str(p)]):
                    with patch.object(w, 'validate_plan'), patch.object(w, 'work') as work:
                        with self.assertRaises(BlockingIOError): w.main()
                        work.assert_not_called()


if __name__ == '__main__': unittest.main()
