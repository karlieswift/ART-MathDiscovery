import sys,unittest,tempfile,pathlib,csv,hashlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'src'))
from n1_benchmark import admit,evaluate,write

class TestAdmission(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.base=pathlib.Path(self.tmp.name)
  self.fasta=self.base/'test.fasta';self.fasta.write_text('>SYNTHETIC_TEST_ACC.1\n'+'ACGT'*100+'\n')
  self.exp=self.base/'historical.tsv';write(self.exp,[],['accession_version','exposure_source'])
  self.row={'source_id':'SYNTHETIC_1','accession_version':'SYNTHETIC_TEST_ACC.1','genome_fasta':str(self.fasta),'genome_sha256':hashlib.sha256(self.fasta.read_bytes()).hexdigest(),'publication_doi':'DOI_TEST_ONLY','array_start0':'20','array_end0':'70','rt_start0':'100','rt_end0':'180','annotation_evidence_uri':'SYNTHETIC_NOT_REAL','prior_contact_status':'CERTIFIED_NEVER_EXPOSED_BEFORE_PRE_FREEZE','family_independence_status':'CERTIFIED_FAMILY_INDEPENDENT','curation_reviewer':'SYNTHETIC_REVIEWER','curation_approval':'APPROVED_BY_HUMAN_CURATOR'}
 def tearDown(self):self.tmp.cleanup()
 def _admit(self,r):
  t=self.base/'input.tsv';write(t,[r],list(r));out=self.base/'decision.tsv';return admit(t,self.exp,out),out
 def test_synthetic_approved(self):
  x,p=self._admit(self.row);self.assertEqual(x['n_admitted_genomes'],1)
 def test_pre_exposure_rejected(self):
  write(self.exp,[{'accession_version':self.row['accession_version'],'exposure_source':'synthetic'}],['accession_version','exposure_source']);x,p=self._admit(self.row);self.assertEqual(x['n_admitted_genomes'],0)
 def test_bad_coords_rejected(self):
  row=dict(self.row);row['array_end0']='3';x,p=self._admit(row);self.assertEqual(x['n_admitted_genomes'],0)
 def test_wrong_checksum_rejected(self):
  row=dict(self.row);row['genome_sha256']='0'*64;x,p=self._admit(row);self.assertEqual(x['n_admitted_genomes'],0)
 def test_evaluation_synthetic_hit(self):
  x,ad=self._admit(self.row)
  truth=self.base/'truth.tsv';write(truth,[{'source_id':self.row['source_id'],'accession_version':self.row['accession_version'],'array_start0':'20','array_end0':'70'}],['source_id','accession_version','array_start0','array_end0'])
  preds=self.base/'preds.tsv';write(preds,[{'accession_version':self.row['accession_version'],'method':'ArmM','start0':'22','end0':'45','rank':'2'},{'accession_version':self.row['accession_version'],'method':'RVT1','start0':'100','end0':'150','rank':'1'}],['accession_version','method','start0','end0','rank'])
  stats=evaluate(ad,truth,preds,self.base/'scores.tsv');self.assertEqual(stats['ArmM']['hits_at_5'],1);self.assertEqual(stats['RVT1']['hits_at_5'],0)
if __name__=='__main__':unittest.main()
