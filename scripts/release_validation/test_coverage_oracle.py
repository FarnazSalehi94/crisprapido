import unittest
from coverage_oracle import enumerate_sites, rc

class CoverageTests(unittest.TestCase):
    def test_exact_strands_coordinates_and_boundaries(self):
        g='ACGTACGTTGCAACGTACGA'
        seq=g+'AGG'+'NNNN'+rc(g+'TGG')
        sites=list(enumerate_sites(seq,g,0,False))
        self.assertEqual([(x['start'],x['end'],x['strand'],x['pam']) for x in sites],[(0,20,'+','AGG'),(30,50,'-','TGG')])
        self.assertEqual(list(enumerate_sites(g+'NGG',g,0,False)),[])
        self.assertEqual(list(enumerate_sites(g+'AG',g,0,False)),[])
    def test_overlapping_and_noncanonical(self):
        sites=list(enumerate_sites('AAAAAGGG','AAAA',0,False))
        self.assertEqual([(s['start'],s['pam']) for s in sites if s['strand']=='+'],[(0,'AGG'),(1,'GGG')])
        self.assertEqual(list(enumerate_sites('ACGT'+'AGC','ACGT',0,False)),[])
    def test_bulges_and_mismatches(self):
        g='ACGTACGTTGCAACGTACGA'
        for target,kind in [(g,'none'),(g[:8]+g[9:],'RNA'),(g[:8]+'T'+g[8:],'DNA')]:
            target=('T' if target[0]!='T' else 'A')+target[1:]
            for seq,strand in [(target+'AGG','+'),(rc(target+'AGG'),'-')]:
                hits=list(enumerate_sites(seq,g))
                self.assertTrue(any(h['strand']==strand and h['bulge']==kind and h['mismatches']==1 for h in hits))
    def test_ambiguous_rejected(self):
        self.assertEqual(list(enumerate_sites('ACNTAGG','ACGT',4,True,0)),[])
if __name__=='__main__': unittest.main()
