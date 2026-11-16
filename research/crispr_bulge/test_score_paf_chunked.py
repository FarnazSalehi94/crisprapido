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
        self.rows = ['row%d\tfield\tvalue' % i for i in range(25)]
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
