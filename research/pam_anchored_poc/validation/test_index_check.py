#!/usr/bin/env python3
import tempfile,unittest,struct
from pathlib import Path
from index_check import check, fasta
class IndexTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();p=Path(self.tmp.name);self.prefix=str(p/'idx');self.ref=p/'ref.fa';self.ref.write_text('>first\nACGT\n>second\nAAN\n')
  Path(self.prefix+'.txt.bin').write_bytes(struct.pack('<I',8)+b'ACGTAAA$');Path(self.prefix+'.pos').write_bytes(struct.pack('<III',0,4,7));Path(self.prefix+'.sna').write_bytes(b''.join(struct.pack('<Q',len(x))+x for x in [b'first',b'second']))
 def tearDown(self):self.tmp.cleanup()
 def test_consistent_index_and_ambiguous_caveat(self):self.assertEqual(check(self.prefix,self.ref)['ambiguous_positions_unverifiable'],1)
 def test_changed_reference_rejected(self):
  self.ref.write_text('>first\nCCGT\n>second\nAAN\n')
  with self.assertRaises(ValueError):check(self.prefix,self.ref)
 def test_reordered_names_rejected(self):
  self.ref.write_text('>other\nACGT\n>second\nAAN\n')
  with self.assertRaises(ValueError):check(self.prefix,self.ref)
 def test_wrapped_fasta_records_match_sequence_and_case(self):
  self.ref.write_text('>first description\naC\ngt\n>second\naa\nn\n')
  self.assertEqual(fasta(self.ref), [['first','ACGT'],['second','AAN']])
  self.assertEqual(check(self.prefix,self.ref)['ambiguous_positions_unverifiable'],1)
 def test_duplicate_and_empty_records_rejected(self):
  for content in ('>same\nA\n>same\nC\n','>first\n>second\nC\n'):
   self.ref.write_text(content)
   with self.assertRaises(ValueError):fasta(self.ref)
 def test_bad_boundary_rejected(self):
  Path(self.prefix+'.pos').write_bytes(struct.pack('<III',0,3,7))
  with self.assertRaises(ValueError):check(self.prefix,self.ref)
if __name__=='__main__':unittest.main()
