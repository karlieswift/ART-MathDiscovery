import sys, unittest, tempfile, json, zipfile, hashlib, random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import discovery_source as s
P=json.loads((s.ROOT/'config/POLICY.json').read_text())
def rec(a='AB123456.1',n=50000,t=123,host=True):
 r={'accession':a,'length':n,'completeness':'COMPLETE','sourceDatabase':'RefSeq',
 'virus':{'taxId':t,'organismName':'Example phage','lineage':[{'taxId':2731619,'name':'Caudoviricetes'}]}}
 if host:r['host']={'taxId':42,'lineage':[{'taxId':2,'name':'Bacteria'}]}
 return r
class Tests(unittest.TestCase):
 def test01_metadata_camel(self):self.assertEqual(s.normalize(rec(),P)[0]['length_bp'],50000)
 def test02_metadata_snake(self):
  r=rec();r['virus']['tax_id']=r['virus'].pop('taxId');r['virus']['organism_name']=r['virus'].pop('organismName');self.assertTrue(s.normalize(r,P)[0])
 def test03_conflicting_alias(self):
  with self.assertRaises(ValueError):s.field({'taxId':1,'tax_id':2},'taxId','tax_id')
 def test04_unknown_host_flag(self):self.assertEqual(s.normalize(rec(host=False),P)[0]['host_status'],'MISSING_HOST_PHAGE_NAME_ONLY')
 def test05_nonbacterial_excluded(self):
  r=rec();r['host']={'taxId':9606,'lineage':[{'taxId':2759}]};self.assertIsNone(s.normalize(r,P)[0])
 def test06_other_virus_excluded(self):
  r=rec();r['virus']['lineage']=[];self.assertIsNone(s.normalize(r,P)[0])
 def test07_incomplete(self):
  r=rec();r['completeness']='PARTIAL';self.assertIsNone(s.normalize(r,P)[0])
 def test08_version_required(self):self.assertIsNone(s.normalize(rec(a='AB123456'),P)[0])
 def test09_lengths(self):
  for n,b in [(20000,'20k_100k'),(99999,'20k_100k'),(100000,'100k_200k'),(200000,'200k_2M'),(2000000,'200k_2M')]:self.assertEqual(s.normalize(rec(n=n),P)[0]['bin'],b)
  self.assertIsNone(s.normalize(rec(n=19999),P)[0]);self.assertIsNone(s.normalize(rec(n=2000001),P)[0])
 def test10_accession_exclude(self):self.assertEqual(s.select([rec()],P,{'AB123456'})[0],[])
 def test11_taxid_representative(self):self.assertEqual(len(s.select([rec(),rec(a='AB123457.1')],P,set())[0]),1)
 def test12_deterministic_input_order(self):
  rr=[rec(a='AB'+str(123400+i)+'.1',t=i+1) for i in range(30)];a=s.select(rr,P,set())[0];rr.reverse();self.assertEqual(a,s.select(rr,P,set())[0])
 def test13_quotas(self):
  p=json.loads(json.dumps(P));p['length_bins'][0]['target']=2
  rr=[rec(a='AB'+str(123400+i)+'.1',t=i+1) for i in range(9)];self.assertEqual(len(s.select(rr,p,set())[0]),2)
 def test14_missing_taxid(self):
  r=rec();r['virus'].pop('taxId');self.assertIsNone(s.normalize(r,P)[0])
 def test15_window_arithmetic(self):
  for n in [0,100,128,2569,5609792]:
   brute=0
   for l in [128,256,512,1024,2048,4096]:
    if n<l:continue
    x=list(range(0,n-l+1,l//4));brute+=2*(len(x)+(x[-1]!=n-l))
   self.assertEqual(s.window_n(n),brute)
 def test16_reverse_complement(self):self.assertEqual(s.reverse_complement('AACGT'),'ACGTT')
 def test17_fasta_no_header(self):
  with self.assertRaises(ValueError):s.parse_fasta('AACGT')
 def test18_wrapped_json_and_jsonl(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'a';p.write_text(json.dumps({'reports':[rec()]}));self.assertEqual(len(list(s.report_objects(p))),1)
   p.write_text(json.dumps(rec())+'\n'+json.dumps(rec())+'\n');self.assertEqual(len(list(s.report_objects(p))),2)
 def test19_zip_md5_and_sequence(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'a.zip';b=b'>AB123456.1\nAACGT\n';fn='ncbi_dataset/data/genomic.fna'
   with zipfile.ZipFile(p,'w') as z:z.writestr(fn,b);z.writestr('md5sum.txt',hashlib.md5(b).hexdigest()+'  '+fn+'\n')
   self.assertEqual(s.zip_genome(p,{'accession':'AB123456.1','length_bp':5}),'AACGT')
 def test20_zip_bad_md5(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'a.zip';fn='ncbi_dataset/data/genomic.fna'
   with zipfile.ZipFile(p,'w') as z:z.writestr(fn,b'>AB123456.1\nAACGT\n');z.writestr('md5sum.txt','0'*32+'  '+fn+'\n')
   with self.assertRaises(ValueError):s.zip_genome(p,{'accession':'AB123456.1','length_bp':5})
 def test21_ambiguity_excluded_without_repair(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'a.zip'
   with zipfile.ZipFile(p,'w') as z:z.writestr('ncbi_dataset/data/genomic.fna','>AB123456.1\nACNT\n')
   with self.assertRaises(s.SourceIneligible):s.zip_genome(p,{'accession':'AB123456.1','length_bp':4})
 def test22_prior240(self):
  rr=s.rows(s.ROOT/'config/EXCLUDE_PRIOR240.tsv');self.assertEqual(len(rr),240);self.assertEqual(sum(int(x['length_bp']) for x in rr),16574805)
  self.assertEqual(len({x['accession'] for x in rr}),240)
 def test23_no_score_or_protein_selection(self):
  r=rec();r['isAnnotated']=True;r['proteinCount']=999999;a=s.select([r],P,set())[0]
  r['isAnnotated']=False;r['proteinCount']=0;self.assertEqual(a,s.select([r],P,set())[0])
 def test24_no_cross_bin_backfill(self):
  rr=[rec(a='AB'+str(123400+i)+'.1',t=i+1) for i in range(260)];self.assertEqual(len(s.select(rr,P,set())[0]),250)
 def test25_invalid_existing_archive_fails(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'a.zip';p.write_text('<html>blocked</html>')
   with self.assertRaises(zipfile.BadZipFile):s.zip_genome(p,{'accession':'AB123456.1','length_bp':5})
if __name__=='__main__':unittest.main(verbosity=2)
