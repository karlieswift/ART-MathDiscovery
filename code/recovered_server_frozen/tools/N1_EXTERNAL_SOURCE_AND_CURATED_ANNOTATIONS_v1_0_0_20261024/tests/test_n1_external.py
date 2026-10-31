import csv,json,pathlib,sys,tempfile,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from curate_public_sources import audit_mars_native, SOURCE_ROWS, FIELDS
from screen_new_phage_metadata import run

class ExternalCurationTests(unittest.TestCase):
 def test_real_mars_native(self):
  with tempfile.TemporaryDirectory() as d:
   r=audit_mars_native(ROOT,pathlib.Path(d))
   self.assertEqual(r['n_verified_anchor_copies'],5)
   self.assertEqual(r['majority_18nt'],'ATATGAATACGTATAATG')
   self.assertFalse(r['primary_benchmark_admissible'])
 def test_source_registry_all_held(self):
  self.assertEqual(len(SOURCE_ROWS),12)
  self.assertTrue(all(x[6].startswith('HOLD') for x in SOURCE_ROWS))
  self.assertTrue(all(x[8].startswith('https://') for x in SOURCE_ROWS))
 def test_new_pool_never_becomes_ART_label(self):
  with tempfile.TemporaryDirectory() as d:
   d=pathlib.Path(d);meta=d/'meta.tsv';exp=d/'exposure.tsv'
   meta.write_text('Accession\tGenome Length (bp)\tDescription\tCheckV quality\nMW248466.1\t266637\tKnown exposed\tHigh-quality\nNEW100003.1\t230000\tNew phage\tHigh-quality\nNEW100004\t250000\tVersion unknown\tComplete\nNEW100005.1\t40000\tSmall\tComplete\n',encoding='utf-8')
   exp.write_text('accession_version\texposure_source\nMW248466.1\tEXPOSED\n',encoding='utf-8')
   r=run(meta,exp,d/'results',min_len=150000)
   self.assertEqual(r['n_priority_queue'],2)
   with (d/'results/N1_UNLABELED_PHAGE_PRIORITY_QUEUE.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
   self.assertEqual(set(r['accession_or_base'] for r in rows),{'NEW100003.1','NEW100004'})
   self.assertTrue(all(r['is_ART_positive']=='UNKNOWN' for r in rows))
   self.assertFalse(next(r for r in rows if r['accession_or_base']=='NEW100004')['exact_version_in_source']=='True')

if __name__=='__main__':unittest.main()
