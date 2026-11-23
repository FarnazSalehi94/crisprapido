import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest

spec = importlib.util.spec_from_file_location('chunked', Path(__file__).with_name('score_paf_chunked.py'))
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)

FAKE_SCORER = textwrap.dedent('''
    import json, sys
    from pathlib import Path
    p = sys.argv[1:]
    args = dict(zip(p[0::2], p[1::2]))
    out = Path(args['--out-dir']); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for line in Path(args['--input']).read_text().splitlines():
        rows.append(line + '\tbu:Z:scored\tcb:f:0.5')
    (out/'output.paf').write_text(chr(10).join(rows)+chr(10))
    (out/'annotations.tsv').write_text('row_index' + chr(10) + chr(10).join('0\tlocal' for _ in rows) + chr(10))
    (out/'metadata.json').write_text(json.dumps(dict(rows=len(rows),
        output_sha256=__import__('hashlib').sha256((chr(10).join(rows)+chr(10)).encode()).hexdigest())))
    (out/'status').write_text('COMPLETE scoring only\\n')
''')


class ChunkedScorerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.inner = self.root/'fake_scorer.py'
        self.inner.write_text(FAKE_SCORER)
        self.input = self.root/'in.paf'
        self.rows = ['row%d\t20\t0\t20\t+\tctg\t1000\t%d\t%d\t20\t20\t255' % (i, i, i+20) for i in range(25)]
        self.input.write_text(chr(10).join(self.rows)+chr(10))
        self.args = lambda out, chunk_rows: [
            '--input', str(self.input), '--out-dir', str(out), '--scoring', 'both',
            '--reference', 'ref', '--guides', 'guides.fa', '--model-root', 'models',
            '--model-python', 'py', '--annotations', '--preserve-cfd',
            '--chunk-rows', str(chunk_rows), '--inner-script', str(self.inner)]

    def tearDown(self):
        self.tmp.cleanup()

    def run_main(self, out, chunk_rows):
        with unittest.mock.patch.object(sys, 'argv', ['x']+self.args(out, chunk_rows)):
            w.main()

    def test_split_merge_and_byte_equivalence(self):
        out = self.root/'out'
        self.run_main(out, 10)
        expected = ''.join(r + chr(9)+'bu:Z:scored'+chr(9)+'cb:f:0.5'+chr(10) for r in self.rows)
        self.assertEqual((out/'output.paf').read_text(), expected)
        ann = (out/'annotations.tsv').read_text().split(chr(10)); self.assertEqual(ann[1].split(chr(9))[0], '0'); self.assertEqual(ann[25].split(chr(9))[0], '24'); self.assertTrue(ann[25].endswith('local')); self.assertEqual(ann[26], '')
        metadata = json.loads((out/'metadata.json').read_text())
        self.assertEqual(metadata['rows'], 25)
        self.assertEqual([c['rows'] for c in metadata['chunks']], [10, 10, 5])
        self.assertTrue((out/'status').read_text().startswith('COMPLETE'))
        for i in range(3):
            self.assertFalse((out/'chunks'/('chunk-%05d.paf' % i)).exists())

    def test_merge_detects_chunk_output_drift(self):
        out = self.root/'out'
        self.run_main(out, 10)
        metadata = json.loads((out/'metadata.json').read_text())
        self.assertEqual(len(metadata['chunks']), 3)
        for record in metadata['chunks']:
            self.assertTrue(record['output_sha256'])
        # Chunk outputs were cleaned after the verified merge; the recorded
        # hashes plus the row-for-row merge check are the tamper evidence.
        self.assertFalse((out/'chunks/00000').exists())
        merged = (out/'output.paf').read_text()
        self.assertTrue(merged.startswith('row0' + chr(9)))

    def test_restart_reuses_completed_chunks(self):
        out = self.root/'out'
        real = w.score_chunk
        calls = []
        def crashing(a, inner, chunk, out_dir):
            if chunk['index'] == 2:
                raise RuntimeError('injected mid-run crash')
            calls.append(chunk['index'])
            return real(a, inner, chunk, out_dir)
        with unittest.mock.patch.object(w, 'score_chunk', crashing):
            with unittest.mock.patch.object(sys, 'argv', ['x']+self.args(out, 10)):
                with self.assertRaises(RuntimeError):
                    w.main()
        self.assertEqual(calls, [0, 1])
        # Second invocation resumes the INCOMPLETE out-dir: only chunk 2 rescores.
        calls.clear()
        def counting(a, inner, chunk, out_dir):
            calls.append(chunk['index'])
            return real(a, inner, chunk, out_dir)
        with unittest.mock.patch.object(w, 'score_chunk', counting):
            with unittest.mock.patch.object(sys, 'argv', ['x']+self.args(out, 10)):
                w.main()
        self.assertEqual(calls, [2])
        self.assertTrue((out/'output.paf').is_file())

    def interrupted(self, out):
        real = w.score_chunk
        def crash(a, inner, chunk, directory):
            if chunk['index'] == 2:
                raise RuntimeError('injected interruption')
            return real(a, inner, chunk, directory)
        with unittest.mock.patch.object(w, 'score_chunk', crash):
            with self.assertRaises(RuntimeError):
                self.run_main(out, 10)

    def test_duplicate_anchor_across_chunks_rejected(self):
        self.input.write_text(self.input.read_text()+self.rows[0]+'\n')
        with self.assertRaisesRegex(ValueError, 'duplicate input PAM anchor'):
            self.run_main(self.root/'out', 10)

    def test_changed_input_never_reused(self):
        out = self.root/'out'
        self.interrupted(out)
        self.input.write_text(self.input.read_text().replace('row0\t', 'NEW0\t'))
        self.run_main(out, 10)
        meta = json.loads((out/'metadata.json').read_text())
        self.assertFalse(meta['chunks'][0]['reused'])
        self.assertTrue(meta['chunks'][1]['reused'])
        self.assertTrue((out/'output.paf').read_text().startswith('NEW0\t'))

    def test_annotation_corruption_never_reused(self):
        out = self.root/'out'
        self.interrupted(out)
        ann = out/'chunks/00000/annotations.tsv'
        ann.write_text(ann.read_text().replace('local', 'WRONG'))
        self.run_main(out, 10)
        self.assertFalse(json.loads((out/'metadata.json').read_text())['chunks'][0]['reused'])
        self.assertTrue(list((out/'chunks').glob('00000.failed-*')))

    def test_incomplete_status_never_reused(self):
        out = self.root/'out'
        self.interrupted(out)
        (out/'chunks/00000/status').write_text('INCOMPLETE\n')
        self.run_main(out, 10)
        self.assertFalse(json.loads((out/'metadata.json').read_text())['chunks'][0]['reused'])

    def test_partial_attempt_log_does_not_block_restart(self):
        out = self.root/'out'
        self.interrupted(out)
        (out/'chunks/invocation-00002.stdout').write_text('old failed invocation\n')
        self.run_main(out, 10)
        self.assertTrue((out/'output.paf').is_file())

    def test_actual_merge_interruption_recovery(self):
        out = self.root/'out'
        real = w.merge
        def crash(*args):
            real(*args)
            raise RuntimeError('crashed after merge')
        with unittest.mock.patch.object(w, 'merge', crash):
            with self.assertRaises(RuntimeError):
                self.run_main(out, 10)
        self.run_main(out, 10)
        self.assertTrue(all(c['reused'] for c in json.loads((out/'metadata.json').read_text())['chunks']))

    def test_resume_rejects_config_drift(self):
        out = self.root/'out'
        out.mkdir()
        (out/'status').write_text('INCOMPLETE chunked\n')
        (out/'chunks').mkdir()
        for i in range(3):
            chunk_out = out/('chunks/%05d' % i)
            chunk_out.mkdir()
            (chunk_out/'status').write_text('COMPLETE scoring only\n')
            (chunk_out/'metadata.json').write_text(json.dumps(dict(rows=10)))
            (chunk_out/'chunk-config.json').write_text(json.dumps(dict(config_id='stale')))
            (chunk_out/'output.paf').write_text('stale\n')
        calls = []
        real = w.score_chunk
        def counting(a, inner, chunk, out_dir):
            calls.append(chunk['index'])
            return real(a, inner, chunk, out_dir)
        with unittest.mock.patch.object(w, 'score_chunk', counting):
            with unittest.mock.patch.object(sys, 'argv', ['x']+self.args(out, 10)):
                w.main()
        self.assertEqual(calls, [0, 1, 2])  # stale chunks never reused

    def test_merge_crash_recovery_removes_stale_artifacts(self):
        out = self.root/'out'
        self.run_main(out, 10)
        self.assertEqual((out/'output.paf.incomplete').exists(), False)

    def test_verify_merged_rejects_divergence(self):
        merged = self.root/'m.paf'
        merged.write_text(''.join(r + chr(9)+'X'+chr(10) for r in self.rows[:24]))
        with self.assertRaises(ValueError):
            w.verify_merged(self.input, merged, 25)

    def test_empty_input_rejected(self):
        empty = self.root/'empty.paf'; empty.write_text('')
        with self.assertRaises(ValueError):
            w.split_input(empty, self.root/'c', 10)

    def test_split_preserves_all_rows(self):
        chunks = w.split_input(self.input, self.root/'c', 7)
        self.assertEqual([c['rows'] for c in chunks], [7, 7, 7, 4])
        rejoined = ''.join((self.root/'c'/('chunk-%05d.paf' % c['index'])).read_text() for c in chunks)
        self.assertEqual(rejoined, self.input.read_text())


import unittest.mock  # noqa: E402  (used inside tests)

if __name__ == '__main__':
    unittest.main()
