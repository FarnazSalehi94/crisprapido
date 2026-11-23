import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('completion', Path(__file__).with_name('scoring_completion.py'))
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)


class CompletionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.r = Path(self.tmp.name)
        self.source = self.r/'in.paf'; self.output = self.r/'out.paf'
        self.ann = self.r/'ann.tsv'; self.spool = self.r/'historical.jsonl'
        self.raw = 'g\t20\t0\t20\t+\tc\t100\t10\t30\t20\t20\t255\tcf:f:0.3000\n'
        self.details = ['0','g','c','10','30','+','A'*20,'A'*20,'TGG','20=','A'*20,'A'*20,'0','0','0','0.3000','available','N/A','scored']
        self.source.write_text(self.raw)
        self.output.write_text(self.raw.rstrip('\n')+'\tcb:f:0.5\tbu:Z:scored\n')
        fields = list(self.details); fields[17] = '0.5'
        self.ann.write_text('\t'.join(w.COLUMNS)+'\n'+'\t'.join(fields)+'\n')
        self.spool.write_text(json.dumps(dict(id=0,raw=self.raw,status='scored',annotation=self.details))+'\n')

    def tearDown(self):
        self.tmp.cleanup()

    def verify(self):
        return w.verify_merged(self.source,self.output,self.ann,self.spool,1,w.digest(self.source),w.digest(self.spool))

    def test_full_verification(self):
        self.assertEqual(self.verify()['availability'], {'scored':1})

    def test_explicit_na(self):
        old = json.loads(self.spool.read_text()); old['status']='unsupported_pam'; old['annotation'][-1]='unsupported_pam'
        self.spool.write_text(json.dumps(old)+'\n')
        self.output.write_text(self.raw.rstrip('\n')+'\tbu:Z:unsupported_pam\n')
        fields = list(self.details); fields[-1]='unsupported_pam'
        self.ann.write_text('\t'.join(w.COLUMNS)+'\n'+'\t'.join(fields)+'\n')
        self.assertEqual(self.verify()['availability'], {'unsupported_pam':1})

    def test_cfd_format_drift_rejected(self):
        self.output.write_text(self.output.read_text().replace('0.3000','0.3'))
        with self.assertRaisesRegex(ValueError,'CFD changed'): self.verify()

    def test_annotation_context_drift_rejected(self):
        self.ann.write_text(self.ann.read_text().replace('TGG','AGG'))
        with self.assertRaisesRegex(ValueError,'annotation/CFD/availability'): self.verify()

    def test_availability_drift_rejected(self):
        self.output.write_text(self.output.read_text().replace('bu:Z:scored','bu:Z:missing'))
        with self.assertRaisesRegex(ValueError,'availability drift'): self.verify()

    def test_nan_rejected(self):
        self.output.write_text(self.output.read_text().replace('cb:f:0.5','cb:f:nan'))
        with self.assertRaisesRegex(ValueError,'invalid CB'): self.verify()

    def test_row_loss_rejected(self):
        self.output.write_text('')
        with self.assertRaisesRegex(ValueError,'row loss'): self.verify()

    def test_extra_rows_rejected(self):
        self.output.write_text(self.output.read_text()*2)
        with self.assertRaisesRegex(ValueError,'extra merged'): self.verify()

    def test_historical_input_drift_rejected(self):
        old=json.loads(self.spool.read_text()); old['raw']=self.raw.replace('cf:f:0.3000','cf:f:0.4')
        self.spool.write_text(json.dumps(old)+'\n')
        with self.assertRaisesRegex(ValueError,'historical conversion order/input'): self.verify()

    def test_historical_digest_rejected(self):
        with self.assertRaisesRegex(ValueError,'historical conversion digest'):
            w.verify_merged(self.source,self.output,self.ann,self.spool,1,w.digest(self.source),'0'*64)

    def test_comparison_tolerance(self):
        a=self.r/'a'; b=self.r/'b'; a.mkdir(); b.mkdir()
        for d in [a,b]:
            (d/'output.paf').write_text(self.output.read_text()); (d/'annotations.tsv').write_text(self.ann.read_text())
        (b/'output.paf').write_text(self.output.read_text().replace('cb:f:0.5','cb:f:0.5000001'))
        (b/'annotations.tsv').write_text(self.ann.read_text().replace('\t0.5\t','\t0.5000001\t'))
        self.assertEqual(w.compare_bundles(a,b)['rows'],1)
        (b/'output.paf').write_text(self.output.read_text().replace('cb:f:0.5','cb:f:0.5001'))
        with self.assertRaisesRegex(ValueError,'tolerance'): w.compare_bundles(a,b)

    def test_no_clobber_publication(self):
        p=self.r/'published.json'; w.publish(p, {'a':1})
        with self.assertRaises(FileExistsError): w.publish(p, {'a':2})
        self.assertEqual(json.loads(p.read_text()), {'a':1})


if __name__ == '__main__': unittest.main()
