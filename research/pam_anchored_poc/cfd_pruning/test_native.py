#!/usr/bin/env python3
"""Small native/CLI regressions. Set the four tool environment variables below."""
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/pam_anchored_poc/cfd_feasibility'))
from study import factor, keep, load_tables, oracle, rc, representatives
from validate import command, metrics, native_rows, paf_rows


class NativeRegressions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cli=os.environ['CRISPR_CFD_CLI']
        cls.native=os.environ['CRISPR_CFD_NATIVE']
        cls.unpruned=os.environ['CRISPR_CFD_UNPRUNED']
        cls.builder=os.environ['COLUMBA_BUILDER']
        cls.tmp=tempfile.TemporaryDirectory(prefix='crispr-cfd-regression-')
        cls.path=Path(cls.tmp.name)
        cls.g='GAAACAGTCGATTTTATCAC'
        cls.tables=load_tables()
        refs={}
        def add(name,target):
            refs[name+'+']=target+'AGG';refs[name+'-']=rc(target+'AGG')
        add('exact',cls.g)
        add('mixed',cls.g[:5]+'C'+cls.g[6:])
        add('alternative',cls.g[:14]+cls.g[15:])
        zero=next((i,t) for i,g in enumerate(cls.g) for t in 'ACGT' if g!=t and factor(g,t,i+1,cls.tables[0])==0)
        i,t=zero
        add('zero',cls.g[:i]+t+cls.g[i+1:])
        cls.panels={}
        for label in ('resolved','N'):
            local=dict(refs)
            if label=='N':
                # Pure D projection leaves a position-1 gap against N. The
                # current scorer treats that gap neutrally BEFORE matrix lookup.
                s=cls.g[0]+'N'+cls.g[1:]+'AGG'
                local['neutralN+']=s;local['neutralN-']=rc(s)
                local['targetN+']='N'+cls.g[1:]+'AGG'
                local['targetN-']=rc(local['targetN+'])
            root=cls.path/label;root.mkdir();ref=root/'ref.fa';guides=root/'guides.fa';idx=root/'index'
            ref.write_text(''.join(f'>{c}\n{s}\n' for c,s in local.items()))
            guides.write_text('>g\n'+cls.g+'\n')
            command([sys.executable,ROOT/'research/pam_anchored_poc/compatibility/index_manifest.py','build',
                     '--builder',cls.builder,'--reference',ref,'--index',idx])
            cls.panels[label]=(local,ref,guides,idx)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_native(self,label,m,z,frac,tau,mode='pruned',pm=None):
        refs,ref,guides,idx=self.panels[label]
        return command([self.native,idx,ref,guides,m,z,frac,'GG',int(z>0),mode,tau,
                        ROOT/'mismatch_scores.txt',pm or ROOT/'pam_scores.txt'])

    def test_native_and_cli_full_precision_equality(self):
        value=.19999999
        pm=self.path/'precision_pam.txt'
        pm.write_text((ROOT/'pam_scores.txt').read_text().replace('GG 1.0','GG '+repr(value)))
        tables=(self.tables[0],dict(self.tables[1],GG=value))
        refs,ref,guides,idx=self.panels['resolved']
        full=representatives(oracle(refs,self.g,'GG',5,0,.75,tables))
        for tau in (0,.1,value,math.nextafter(value,-math.inf),math.nextafter(value,math.inf),.2,.5):
            expected={k:v for k,v in full.items() if keep(v,tau)}
            a=self.run_native('resolved',5,0,.75,tau,'unpruned',pm)
            b=self.run_native('resolved',5,0,.75,tau,'pruned',pm)
            native_rows(a.stdout,expected);native_rows(b.stdout,expected)
            out=command([self.cli,'--search-backend','pam-anchored','--pam-search-bin',self.unpruned,
                         '--columba-index',idx,'-r',ref,'-G',guides,'-t','1','-b','0','--pam-scores',pm,
                         '--min-cfd',repr(tau)])
            paf_rows(out.stdout,expected)
            if tau==.2:
                self.assertEqual(len(expected),0)
                self.assertGreater(metrics(b.stderr)['pruned'],0)

    def test_mixed_only_unknown_representatives_survive(self):
        refs,_,_,_=self.panels['resolved']
        full=representatives(oracle(refs,self.g,'GG',0,2,.75,self.tables))
        unknown={k:v for k,v in full.items() if v.score is None}
        self.assertTrue(any(k[0].startswith('mixed') for k in unknown))
        for tau in (0,.1,.2,.5,1):
            r=self.run_native('resolved',0,2,.75,tau)
            native_rows(r.stdout,{k:v for k,v in full.items() if keep(v,tau)})
            self.assertGreaterEqual(metrics(r.stderr)['unavailable'],len(unknown))

    def test_neutral_position_one_gap_with_n(self):
        refs,_,_,_=self.panels['N']
        full=representatives(oracle(refs,self.g,'GG',5,2,.75,self.tables))
        neutral=[x for k,x in full.items() if k[0].startswith('neutralN')]
        self.assertEqual(len(neutral),2)
        self.assertTrue(all(x.score==1 for x in neutral))
        r=self.run_native('N',5,2,.75,.5)
        native_rows(r.stdout,{k:v for k,v in full.items() if keep(v,.5)})
        self.assertEqual(metrics(r.stderr)['pruned'],0)
        self.assertFalse(metrics(r.stderr)['bound_eligible'])

    def test_research_binary_is_not_integrated_backend_protocol(self):
        _,ref,guides,idx=self.panels['resolved']
        r=subprocess.run(list(map(str,[self.cli,'--search-backend','pam-anchored',
            '--pam-search-bin',self.native,'--columba-index',idx,'-r',ref,'-G',guides,
            '-t','1','--min-cfd','.1'])),cwd=ROOT,capture_output=True,text=True,timeout=30)
        self.assertNotEqual(r.returncode,0)
        self.assertEqual(r.stdout,'')
        self.assertIn('experimental backend failed',r.stderr)

    def test_zero_scores_are_not_unavailable(self):
        refs,_,_,_=self.panels['resolved']
        full=representatives(oracle(refs,self.g,'GG',5,0,.75,self.tables))
        zeros=[v for k,v in full.items() if k[0].startswith('zero')]
        self.assertTrue(zeros)
        self.assertTrue(all(v.score==0 for v in zeros))
        for tau in (0,.1,.2,.5):
            r=self.run_native('resolved',5,0,.75,tau)
            native_rows(r.stdout,{k:v for k,v in full.items() if keep(v,tau)})


if __name__=='__main__':
    unittest.main()
