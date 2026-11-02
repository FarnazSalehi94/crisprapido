import json,subprocess,tempfile,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent

def sam(pos,cigar='20M',flag=0):return f'g\t{flag}\tchr\t{pos+1}\t60\t{cigar}\t*\t0\t0\t*\t*\tNM:i:0\tAS:i:0\n'
def paf(start,end,strand='+'):return f'g\t20\t0\t20\t{strand}\tchr\t1000\t{start}\t{end}\t20\t20\t255\n'
class CandidateComparisonTests(unittest.TestCase):
 def test_exact_key_and_start_are_distinguished(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);(p/'old.sam').write_text(sam(10)+sam(40)+sam(70));(p/'new.sam').write_text(sam(10)+sam(40,'19M1I')+sam(100))
   subprocess.run(['python3',str(HERE/'compare_columba_candidate_sets.py'),'--old',str(p/'old.sam'),'--new',str(p/'new.sam'),'--output',str(p/'out.json')],check=True,stdout=subprocess.DEVNULL)
   x=json.loads((p/'out.json').read_text());self.assertEqual((x['old_keys_not_exactly_present'],x['old_starts_not_present_in_new']),(2,1))
 def test_strand_and_actual_span_control_overlap(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);(p/'old.paf').write_text(paf(20,40)+paf(100,120,'-')+paf(200,220))
   (p/'new.sam').write_text(sam(39)+sam(110,flag=16)+sam(180,'20M'))
   subprocess.run(['python3',str(HERE/'compare_prior_anchor_overlap.py'),'--old-paf',str(p/'old.paf'),'--new-sam',str(p/'new.sam'),'--output',str(p/'out.json')],check=True,stdout=subprocess.DEVNULL)
   x=json.loads((p/'out.json').read_text());self.assertEqual((x['old_rows_with_overlapping_candidate'],x['old_rows_without_overlapping_candidate']),(2,1))
if __name__=='__main__':unittest.main()
