import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest

s=importlib.util.spec_from_file_location('publication',Path(__file__).with_name('publish_scoring_results.py'))
w=importlib.util.module_from_spec(s);s.loader.exec_module(w)


class PublicationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.r=Path(self.tmp.name)
        subprocess.run(['git','init','-q','-b','task',str(self.r)],check=True)
        w.git(self.r,'config','user.name','Test');w.git(self.r,'config','user.email','test@example.invalid')
        (self.r/'docs').mkdir();(self.r/w.HANDOFF).write_text(w.HEADER+'\nBase history.\n')
        (self.r/'other').write_text('base\n')
        w.git(self.r,'add','.');w.git(self.r,'commit','-qm','base')
        w.git(self.r,'remote','add','origin','https://example.invalid/project.git')
        self.head=w.git(self.r,'rev-parse','HEAD')
        (self.r/w.HANDOFF).write_text(w.HEADER+'\nUSER UNCOMMITTED NOTE\n\nBase history.\n')
        (self.r/'other').write_text('USER OTHER CHANGE\n')

    def tearDown(self):self.tmp.cleanup()

    def commit(self):
        return w.commit_results(self.r,self.head,'task','https://example.invalid/project.git','actual result\n','\nOUR RESULT\n')

    def test_preserves_unrelated_changes(self):
        commit=self.commit()
        self.assertNotEqual(commit,self.head)
        committed=w.git(self.r,'show','HEAD:'+w.HANDOFF)
        self.assertIn('OUR RESULT',committed);self.assertNotIn('USER UNCOMMITTED',committed)
        self.assertIn('USER UNCOMMITTED',(self.r/w.HANDOFF).read_text())
        self.assertEqual((self.r/'other').read_text(),'USER OTHER CHANGE\n')
        self.assertEqual(w.git(self.r,'diff','--cached','--name-only'),'')
        self.assertEqual(set(w.git(self.r,'diff','--name-only').splitlines()),{w.HANDOFF,'other'})
        self.assertEqual(w.git(self.r,'show','HEAD:'+w.RESULT),'actual result')

    def test_branch_drift_fails_before_writes(self):
        w.git(self.r,'switch','-qc','other-branch')
        with self.assertRaisesRegex(ValueError,'branch changed'):self.commit()
        self.assertFalse((self.r/w.RESULT).exists())

    def test_head_drift_fails_before_writes(self):
        w.git(self.r,'commit','--allow-empty','-qm','other commit')
        with self.assertRaisesRegex(ValueError,'HEAD changed'):self.commit()
        self.assertFalse((self.r/w.RESULT).exists())

    def test_staged_change_fails_before_writes(self):
        w.git(self.r,'add','other')
        with self.assertRaisesRegex(ValueError,'staged changes'):self.commit()
        self.assertFalse((self.r/w.RESULT).exists())

    def test_existing_result_not_overwritten(self):
        (self.r/w.RESULT).write_text('USER RESULT\n')
        with self.assertRaisesRegex(ValueError,'already exists'):self.commit()
        self.assertEqual((self.r/w.RESULT).read_text(),'USER RESULT\n')

    def test_remote_drift_fails_before_writes(self):
        w.git(self.r,'remote','set-url','origin','https://example.invalid/other.git')
        with self.assertRaisesRegex(ValueError,'remote changed'):self.commit()
        self.assertFalse((self.r/w.RESULT).exists())


if __name__=='__main__':unittest.main()
