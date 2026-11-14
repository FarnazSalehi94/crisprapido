import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).parent))
import review_quarantine as deletion
import validate_replacement as gate


class RecoveryGateTests(unittest.TestCase):
    def test_apply_bad_review_digest_preserves_files(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'review.json';p.write_text('{}')
            with self.assertRaisesRegex(ValueError,'manifest changed'):
                deletion.apply(p,'0'*64)
            self.assertTrue(p.exists());self.assertFalse(p.with_suffix('.deletion.json').exists())

    def test_apply_refuses_changed_live_review_before_unlink(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'review.json';m=dict(root=d,assembly='A',validation=d,monitor=d)
            p.write_text(json.dumps(m))
            with patch.object(deletion,'review',return_value=dict(m,changed=True)):
                with self.assertRaisesRegex(ValueError,'live evidence differs'):
                    deletion.apply(p,deletion.digest(p))
            self.assertFalse(p.with_suffix('.deletion.json').exists())

    def test_apply_regular_files_only_logs_removed_bytes(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);q=root/'q';q.mkdir();f=q/'A.bwt';f.write_bytes(b'test')
            st=f.stat()
            m=dict(root=d,assembly='A',validation=d,monitor=d,quarantine=str(q),
                   entries=[dict(name=f.name,inode=st.st_ino,bytes=4,mtime_ns=st.st_mtime_ns)])
            p=root/'review.json';p.write_text(json.dumps(m))
            with patch.object(deletion,'review',return_value=m):
                result=deletion.apply(p,deletion.digest(p))
            self.assertEqual(result['status'],'COMPLETE');self.assertEqual(result['logical_removed_bytes'],4)
            self.assertFalse(q.exists());self.assertTrue(p.exists())

    def test_symlink_swap_is_detected_not_followed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);q=root/'q';q.mkdir();target=root/'healthy.fa';target.write_text('ACGT')
            f=q/'A.bwt';f.symlink_to(target);st=f.lstat()
            m=dict(root=d,assembly='A',validation=d,monitor=d,quarantine=str(q),
                   entries=[dict(name=f.name,inode=st.st_ino,bytes=st.st_size,mtime_ns=st.st_mtime_ns)])
            p=root/'review.json';p.write_text(json.dumps(m))
            with patch.object(deletion,'review',return_value=m):
                with self.assertRaisesRegex(ValueError,'quarantine changed'):
                    deletion.apply(p,deletion.digest(p))
            self.assertEqual(target.read_text(),'ACGT')
            self.assertEqual(json.loads(p.with_suffix('.deletion.json').read_text())['status'],'PARTIAL')

    def test_validation_rejects_same_build_and_recheck_host(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);index=root/'indices/A';index.mkdir(parents=True)
            manifest=dict(assembly='A',status='TRUSTED_CROSS_NODE',
                          cross_node_recheck=dict(result='PASS',all_components_match=True,host='n1',build_host='n1'))
            (index/'completion.json').write_text(json.dumps(manifest));(index/'build-host.txt').write_text('BUILD_HOST=n1\n')
            with self.assertRaisesRegex(ValueError,'cross-host'):
                gate.validate(root,'A',root/'comparator','0'*64,root/'out.json')
            self.assertFalse((root/'out.json').exists())


if __name__=='__main__':unittest.main()
