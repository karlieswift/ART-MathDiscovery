import sys,unittest,tempfile,pathlib,csv
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'src'))
from n3_analysis import calculate
from n4_analysis import analyze
ROOT=pathlib.Path(__file__).resolve().parents[1]
class TestRecalculation(unittest.TestCase):
 def test_n3_real_archived_data(self):
  with tempfile.TemporaryDirectory() as d:
   r=calculate(ROOT/'data/N3_REGION_RANKS_6_NEW_MATRICES.tsv.gz',d)
   self.assertEqual(r['n_combinations'],18);self.assertTrue(r['all_identical_at_frozen_threshold'])
 def test_n4_real_archived_pairing(self):
  with tempfile.TemporaryDirectory() as d:
   r=analyze(ROOT/'data/N4_Review69_SANITIZED.tsv',ROOT/'data/N4_Pairs2144_SANITIZED.tsv',d)
   self.assertEqual(r['n_pairs'],2144);self.assertEqual(r['strict_positive_count'],24)
if __name__=='__main__':unittest.main()
