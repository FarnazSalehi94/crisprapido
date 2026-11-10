"""Opt-in CLI regressions; point CRISPR_SCORING_CLI at a development build."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

CLI = os.environ.get('CRISPR_SCORING_CLI')
MODEL_PYTHON = os.environ.get('CRISPR_MODEL_PYTHON')
MODEL_ROOT = os.environ.get('CRISPR_MODEL_ROOT')
GUIDE = 'GAAACAGTCGATTTTATCAC'


@unittest.skipUnless(CLI, 'requires isolated CRISPR_SCORING_CLI development binary')
class CliTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.ref = self.root/'ref.fa'
        self.ref.write_text('>chr\n'+GUIDE+'AGG\n')
        self.paf = self.root/'input.paf'
        self.paf.write_text('Guide\t20\t0\t20\t+\tchr\t23\t0\t20\t20\t20\t255\tnm:i:0\tng:i:0\tbs:i:0\tcg:Z:20=\tcf:f:1.0000\n')
        self.dummy = self.root/'models'; self.dummy.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def run_cli(self, opts, saved=True, success=True):
        command = [CLI, '-r', str(self.ref), '-g', GUIDE]
        if saved:
            command += ['--score-paf', str(self.paf)]
        result = subprocess.run(command+opts, capture_output=True, timeout=60)
        self.assertEqual(result.returncode == 0, success, result.stderr.decode())
        return result

    def neural(self, mode, root, real=False):
        return ['--scoring', mode, '--scoring-output-dir', str(root),
                '--crispr-bulge-python', MODEL_PYTHON if real else str(Path(sys.executable).absolute()),
                '--crispr-bulge-model-root', MODEL_ROOT if real else str(self.dummy)]

    def test_cfd_default_pass_through_and_option_errors(self):
        for options in ([], ['--scoring', 'cfd']):
            self.assertEqual(self.run_cli(options).stdout, self.paf.read_bytes())
        for options in (['--scoring', 'invalid'], ['--scoring', 'both'], ['--min-cfd', '.5'],
                        ['--scoring-batch-size', '0'], ['--scoring-worker-timeout', '0']):
            self.assertFalse(self.run_cli(options, success=False).stdout)

    def test_unsupported_is_na_and_existing_cfd_or_missing_cfd_preserved(self):
        cases = [('DNA2', GUIDE[:5]+'AA'+GUIDE[5:], '5=2D15=', 20, 22,
                  'alignment_exceeds_24_columns', GUIDE[:5]+'--'+GUIDE[5:]),
                 ('RNA2', GUIDE[:5]+GUIDE[7:], '5=2I13=', 18, 20,
                  'unsupported_multibase_bulge', GUIDE),
                 ('mixed', GUIDE[:5]+'A'+GUIDE[6:], '5=1I1D14=', 19, 21,
                  'unsupported_mixed_bulge', GUIDE[:6]+'-'+GUIDE[6:])]
        for name, target, cigar, matches, block, reason, aligned in cases:
            self.ref.write_text('>chr\n'+target+'AGG\n')
            raw = f'Guide\t20\t0\t20\t+\tchr\t{len(target)+3}\t0\t{len(target)}\t{matches}\t{block}\t255\tnm:i:0\tng:i:1\tbs:i:2\tcg:Z:{cigar}'
            for mode in ('crispr-bulge', 'both'):
                for present in (False, True):
                    self.paf.write_text(raw+('\tcf:f:0.0000' if present else '')+'\n')
                    root = self.root/f'run-{name}-{mode}-{present}'
                    result = self.run_cli(self.neural(mode, root))
                    self.assertTrue(result.stdout.startswith(self.paf.read_bytes().rstrip()+b'\t'))
                    self.assertNotIn(b'cb:f:', result.stdout)
                    with (root/'scoring/annotations.tsv').open() as f:
                        annotation = list(csv.DictReader(f, delimiter='\t'))
                    self.assertEqual(len(annotation), 1)
                    self.assertEqual(annotation[0]['crispr_bulge'], 'N/A')
                    self.assertEqual(annotation[0]['crispr_bulge_reason'], reason)
                    self.assertEqual(annotation[0]['cfd'], '0.0000' if present else 'N/A')
                    self.assertEqual(annotation[0]['aligned_guide'], aligned)
                    self.assertIsNone(json.loads((root/'scoring/metadata.json').read_text())['worker'])

    def test_empty_input_is_complete_without_loading_models(self):
        self.paf.write_text('')
        root = self.root/'empty'
        self.assertEqual(self.run_cli(self.neural('both', root)).stdout, b'')
        self.assertEqual(len((root/'scoring/annotations.tsv').read_text().splitlines()), 1)
        self.assertIsNone(json.loads((root/'scoring/metadata.json').read_text())['worker'])

    def test_worker_failure_and_existing_directory_publish_no_stdout(self):
        root = self.root/'failed'
        result = self.run_cli(self.neural('both', root), success=False)
        self.assertFalse(result.stdout)
        self.assertEqual((root/'status').read_text(), 'INCOMPLETE\n')
        self.assertFalse((root/'scoring/output.paf').exists())
        self.assertFalse(self.run_cli(self.neural('both', root), success=False).stdout)

    @unittest.skipUnless(MODEL_PYTHON and MODEL_ROOT, 'requires optional pinned model environment')
    def test_live_search_modes_batch_reuse_and_default_compatibility(self):
        sam = self.root/'input.sam'
        sam.write_text('guide\t0\tchr\t1\t255\t20M\t*\t0\t0\t*\t*\tNM:i:0\tAS:i:0\n')
        search = ['--columba-sam', str(sam), '-t', '1', '-m', '0', '-b', '0']
        default = self.run_cli(search, saved=False).stdout
        self.assertTrue(default)
        self.assertEqual(default, self.run_cli(search+['--scoring', 'cfd'], saved=False).stdout)
        old = os.environ.get('CRISPR_OLD_CLI')
        if old:
            expected = subprocess.run([old, '-r', str(self.ref), '-g', GUIDE]+search,
                                      capture_output=True, check=True).stdout
            self.assertEqual(default, expected)
        for mode in ('crispr-bulge', 'both'):
            root = self.root/mode
            result = self.run_cli(search+self.neural(mode, root, real=True), saved=False)
            self.assertTrue(result.stdout.startswith(default.rstrip()+b'\tcb:f:'))
            self.assertEqual(json.loads((root/'scoring/worker-metrics.json').read_text())['model_load_count'], 5)
            self.assertEqual((root/'input.paf').read_bytes(), default)
        # Multi-guide reporting is spooled as one model invocation, not per guide.
        guides = self.root/'guides.fa'
        guides.write_text('>one\n'+GUIDE+'\n>two\n'+GUIDE+'\n')
        sam.write_text(''.join(f'{name}\t0\tchr\t1\t255\t20M\t*\t0\t0\t*\t*\tNM:i:0\tAS:i:0\n'
                               for name in ('one', 'two')))
        root = self.root/'batch'
        result = subprocess.run([CLI, '-r', str(self.ref), '-G', str(guides)]+search+
                                self.neural('both', root, real=True), capture_output=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        self.assertEqual([s.split(b'\t')[0] for s in result.stdout.splitlines()], [b'one', b'two'])
        worker = json.loads((root/'scoring/worker-metrics.json').read_text())
        self.assertEqual(worker['model_load_count'], 5)
        self.assertEqual(sum(b['rows'] for b in worker['batches']), 2)
        # A timed-out worker is a failed run, never thousands of N/A predictions.
        root = self.root/'timeout'
        failed = self.run_cli(self.neural('both', root, real=True)+['--scoring-worker-timeout', '1'], success=False)
        self.assertFalse(failed.stdout)
        self.assertEqual((root/'status').read_text(), 'INCOMPLETE\n')
        self.assertFalse((root/'scoring/output.paf').exists())


if __name__ == '__main__':
    unittest.main()
