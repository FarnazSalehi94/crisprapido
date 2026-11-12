#!/usr/bin/env python3
"""Fast stdlib-only tests; TensorFlow is deliberately not required."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from adapter import IndexedFasta, convert, encode_tokens, rc

GUIDE = 'GAAACAGTCGATTTTATCAC'
SCRIPT = Path(__file__).with_name('score_paf.py')


def paf(target, cigar, strand='+', name='g', contig='ref', pam='AGG', guide=GUIDE):
    seq = target+pam
    seq = seq if strand == '+' else rc(seq)
    lo, hi = (0,len(target)) if strand == '+' else (len(pam),len(seq))
    return seq, '\t'.join(map(str,[name,len(guide),0,len(guide),strand,contig,len(seq),lo,hi,0,20,255]))+'\tcg:Z:'+cigar+'\tcf:f:0.1234\n'


def changed(seq):
    return ''.join('C' if x != 'C' else 'A' for x in seq)


class AdapterTests(unittest.TestCase):
    def converted(self,target,cigar,strand='+',pam='AGG',guide=GUIDE):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'ref.fa'
            seq,row=paf(target,cigar,strand,pam=pam,guide=guide)
            path.write_text('>ref\n'+seq+'\n')
            ref=IndexedFasta(path)
            try:return convert(row,{'g':guide},ref)
            finally:ref.close()

    def test_supported_inputs_on_both_strands(self):
        cases=[(GUIDE,'20='),(GUIDE[:5]+'T'+GUIDE[5:],'5=1D15='),
               (GUIDE[:5]+GUIDE[6:],'5=1I14='),
               (changed(GUIDE[:4])+'T'+GUIDE[4:],'4X1D16='),
               (changed(GUIDE[:6])+GUIDE[6:],'6X14=')]
        for t,cg in cases:
            for s in '+-':
                _,status,features=self.converted(t,cg,s)
                self.assertEqual(status,'scored')
                self.assertEqual(len(encode_tokens([features])[0]),24)
                self.assertTrue(features['guide'].endswith('NGG'))
                self.assertTrue(features['target'].endswith('AGG'))
                self.assertEqual(features['guide'].replace('-','')[:-3],GUIDE)

    def test_nrr_all_mismatch_counts_gap_types_and_strands(self):
        for suffix in ('GG', 'AG', 'GA', 'AA'):
            for first in 'ACGT':
                pam = first+suffix
                for m in range(7):
                    prefix = (str(m)+'X') if m else ''
                    for gap in ('none', 'DNA', 'RNA'):
                        if gap == 'none':
                            target = changed(GUIDE[:m])+GUIDE[m:]
                            cigar = prefix+str(20-m)+'='
                        elif gap == 'DNA':
                            target = changed(GUIDE[:m])+'T'+GUIDE[m:]
                            cigar = prefix+'1D'+str(20-m)+'='
                        else:
                            target = changed(GUIDE[:m])+GUIDE[m+1:]
                            cigar = prefix+'1I'+str(19-m)+'='
                        for strand in '+-':
                            _, status, features = self.converted(target,cigar,strand,pam=pam)
                            expected = 'scored' if gap == 'none' or m <= 4 else 'outside_documented_mismatch_domain'
                            self.assertEqual(status, expected)
                            if features:
                                self.assertTrue(features['guide'].endswith('NGG'))
                                self.assertTrue(features['target'].endswith(pam))
                                tokens = encode_tokens([features])[0]
                                self.assertEqual(tokens[-2:], [10+'ACGT'.index(x) for x in suffix])
                for target,cigar,reason in [(GUIDE[:5]+'TT'+GUIDE[5:],'5=2D15=','alignment_exceeds_24_columns'), (GUIDE[:5]+GUIDE[7:],'5=2I13=','unsupported_multibase_bulge')]:
                    self.assertEqual(self.converted(target,cigar,pam=pam)[1],reason)

    def test_unsupported_domain_is_not_fabricated(self):
        cases=[(GUIDE[:5]+'TT'+GUIDE[5:],'5=2D15=','alignment_exceeds_24_columns'),
               (GUIDE[:5]+GUIDE[7:],'5=2I13=','unsupported_multibase_bulge'),
               (GUIDE[:5]+'T'+GUIDE[6:],'5=1I1D14=','unsupported_mixed_bulge'),
               (changed(GUIDE[:5])+'T'+GUIDE[5:],'5X1D15=','outside_documented_mismatch_domain'),
               (changed(GUIDE[:7])+GUIDE[7:],'7X13=','outside_documented_mismatch_domain'),
               ('N'+GUIDE[1:],'1X19=','ambiguous_reference')]
        for target,cigar,expected in cases:
            for strand in '+-':
                _,status,features=self.converted(target,cigar,strand)
                self.assertEqual(status,expected)
                self.assertIsNone(features)
        self.assertEqual(self.converted(GUIDE,'20=',pam='NGG')[1],'ambiguous_pam')
        self.assertEqual(self.converted(GUIDE,'20=',pam='AAC')[1],'unsupported_pam')
        self.assertEqual(self.converted(GUIDE,'20=',guide='N'+GUIDE[1:])[1],'ambiguous_guide')

    def test_independent_encoding_padding_and_structural_n(self):
        tokens=encode_tokens([{'guide':GUIDE+'NGG','target':GUIDE+'AGG'}])[0]
        self.assertEqual(tokens[0],24)  # upstream GAP padding, not zero/AA padding
        self.assertEqual(tokens[-3:], [0,12,12])  # N:A -> A:A, G:G, G:G
        with self.assertRaises(ValueError):
            encode_tokens([{'guide':GUIDE[:5]+'--'+GUIDE[5:]+'NGG','target':GUIDE[:5]+'TT'+GUIDE[5:]+'AGG'}])
        with self.assertRaises(ValueError):
            encode_tokens([{'guide':GUIDE+'NGG','target':'N'+GUIDE[1:]+'AGG'}])

    def test_indexed_fasta_wrapping_crlf_and_last_line(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'ref.fa'
            p.write_bytes(b'>a\r\nACGT\r\nNNAC\r\nGT\r\n>b\nTACG')
            r=IndexedFasta(p)
            try:
                self.assertEqual(r.fetch('a',2,10),'GTNNACGT')
                self.assertEqual(r.fetch('b',0,4),'TACG')
                with self.assertRaises(ValueError):r.fetch('a',0,11)
            finally:r.close()
            p.write_text('>a\nACGT\nAC\nGT\n')
            with self.assertRaises(ValueError):IndexedFasta(p)

    def test_cigar_assertions_and_consumption_fail_clearly(self):
        with self.assertRaises(ValueError):self.converted('C'+GUIDE[1:],'20=')
        with self.assertRaises(ValueError):self.converted(GUIDE,'19=')
        with self.assertRaises(ValueError):self.converted(GUIDE,'20=Q')

    def test_default_cfd_mode_is_byte_identical_without_tensorflow(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);seq,row=paf(GUIDE,'20=');p=d/'input.paf'
            p.write_bytes(row.replace('\n','\r\n').encode())
            r=subprocess.run([sys.executable,SCRIPT,'--input',p,'--out-dir',d/'out'],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual(p.read_bytes(),(d/'out/output.paf').read_bytes())
            self.assertIsNone(json.loads((d/'out/metadata.json').read_text())['worker'])
            again=subprocess.run([sys.executable,SCRIPT,'--input',p,'--out-dir',d/'out'],capture_output=True)
            self.assertNotEqual(again.returncode,0)

    def test_model_checksum_gate_rejects_changed_weights(self):
        from worker import verify_models
        from adapter import digest
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);weight=root/'model.h5';weight.write_bytes(b'checksum-fixture-only')
            spec={'models':[{'path':'model.h5','sha256':digest(weight)}]}
            self.assertEqual(verify_models(root,spec),[weight])
            weight.write_bytes(b'changed')
            with self.assertRaises(ValueError):verify_models(root,spec)

    def test_duplicate_or_wrong_reference_leaves_incomplete_result(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);seq,row=paf(GUIDE,'20=');ref=d/'ref.fa';ref.write_text('>ref\n'+seq+'\n')
            guides=d/'guides.fa';guides.write_text('>g\n'+GUIDE+'\n');p=d/'input.paf';p.write_text(row+row)
            common=[sys.executable,SCRIPT,'--input',p,'--scoring','both','--reference',ref,
                    '--guides',guides,'--model-root',d,'--model-python',sys.executable]
            for name,extra in [('duplicate',[]),('wrong_reference',['--expected-reference-sha256','0'*64])]:
                dest=d/name
                r=subprocess.run(common+['--out-dir',dest]+extra,capture_output=True,text=True)
                self.assertNotEqual(r.returncode,0)
                self.assertFalse((dest/'output.paf').exists())
                self.assertEqual((dest/'status').read_text(),'INCOMPLETE\n')

    def test_worker_failure_is_not_published_as_unavailable_predictions(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);seq,row=paf(GUIDE,'20=');ref=d/'ref.fa';ref.write_text('>ref\n'+seq+'\n')
            guides=d/'guides.fa';guides.write_text('>g\n'+GUIDE+'\n');p=d/'input.paf';p.write_text(row)
            failure=d/'failing-interpreter';failure.write_text('#!/bin/sh\nexit 7\n');failure.chmod(0o700)
            r=subprocess.run([sys.executable,SCRIPT,'--input',p,'--out-dir',d/'out','--scoring','both',
                '--reference',ref,'--guides',guides,'--model-root',d,'--model-python',failure],capture_output=True,text=True)
            self.assertNotEqual(r.returncode,0)
            self.assertIn('worker failed',r.stderr)
            self.assertFalse((d/'out/output.paf').exists())
            self.assertEqual((d/'out/status').read_text(),'INCOMPLETE\n')

    def test_unsupported_rows_preserved_and_worker_not_loaded(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);target=GUIDE[:5]+'TT'+GUIDE[5:];seq,row=paf(target,'5=2D15=')
            ref=d/'ref.fa';ref.write_text('>ref\n'+seq+'\n')
            p=d/'input.paf';p.write_text(row)
            for mode in ('both','crispr-bulge'):
                result=d/mode
                r=subprocess.run([sys.executable,SCRIPT,'--input',p,'--out-dir',result,'--scoring',mode,
                    '--reference',ref,'--guide',GUIDE,'--model-root',d,'--model-python',sys.executable,
                    '--guides',d/'not_used'],capture_output=True,text=True)
                # Mutually exclusive guides/guide must be rejected before any model launch.
                self.assertNotEqual(r.returncode,0)
                guidefile=d/'guides.fa';guidefile.write_text('>g\n'+GUIDE+'\n')
                r=subprocess.run([sys.executable,SCRIPT,'--input',p,'--out-dir',result,'--scoring',mode,
                    '--reference',ref,'--guides',guidefile,'--model-root',d,'--model-python',sys.executable],capture_output=True,text=True)
                self.assertEqual(r.returncode,0,r.stderr)
                text=(result/'output.paf').read_text()
                self.assertIn('bu:Z:alignment_exceeds_24_columns',text)
                self.assertNotIn('cb:f:',text)
                self.assertEqual('cf:f:0.1234' in text,mode=='both')
                self.assertIsNone(json.loads((result/'metadata.json').read_text())['worker'])


if __name__=='__main__':
    unittest.main()
