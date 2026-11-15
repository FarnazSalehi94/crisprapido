import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('prepare', Path(__file__).with_name('prepare_recovery_wave.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SelectionTests(unittest.TestCase):
    def test_disjoint_and_never_rebuild_existing_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            rows = [dict(assembly='A'+str(i), status='INDEX_EIO') for i in range(466)]
            rows[0]['status'] = 'READY_BOUND'
            binding = root/'runs/sbds_ngg_nrr_v1/00_manifest/bindings/A1/binding.json'
            binding.parent.mkdir(parents=True); binding.touch()
            (root/'indices_rebuilt_v1/A2').mkdir(parents=True)
            promoted = root/'indices/A3/promotion.json'
            promoted.parent.mkdir(parents=True); promoted.touch()
            self.assertEqual(module.choose(rows, root, {'A4'}, 2), ['A5', 'A6'])

    def test_duplicate_cohort_rejected(self):
        with self.assertRaises(ValueError):
            module.choose([dict(assembly='A', status='INDEX_EIO')]*466, Path('/unused'), set(), 1)


if __name__ == '__main__':
    unittest.main()
