#!/usr/bin/env python3
"""Executable guard tests. Synthetic fixtures are never scientific evidence."""
import argparse,contextlib,importlib.metadata,io,json,sys,tempfile,types,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
import j2_pipeline as j

class Guards(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory(prefix='j2_guard_');self.root=Path(self.temp.name)
 def tearDown(self):self.temp.cleanup()
 def make_run(self,failed=False):
  r=self.root/'fixture';f=r/'frozen';s=r/'sealed';o=r/'scores';f.mkdir(parents=True);s.mkdir();o.mkdir()
  (f/'driver').mkdir();(f/'driver/j2_pipeline.py').write_bytes(Path(j.__file__).read_bytes())
  pol=dict(j.POLICY,analysis_kind='TEST_FIXTURE_ONLY_NOT_SCIENTIFIC',package_versions={x:importlib.metadata.version(x) for x in ['numpy','pandas']})
  j.write_json(f/'POLICY.json',pol)
  j.write_json(f/'TOOLCHAIN.json',{'python':{'path':sys.executable,'sha256':j.sha(Path(sys.executable).resolve())}})
  j.write_tsv(f/'SCORING_MANIFEST.tsv',[{'opaque_id':'J2_00001','length_bp':10000}])
  j.write_tsv(f/'FROZEN_FILE_SHA256.tsv',j.hashes(r,list(f.rglob('*'))))
  j.write_tsv(s/'IDENTITY_KEY.tsv',[{'opaque_id':'J2_00001','sample_id':'fixture'}])
  j.write_tsv(s/'TRUTH.tsv',[{'sample_id':'fixture','target_id':'fixture_unit_1','label_class':'positive','truth_start0':100,'truth_end0':200,'array_start0':100,'array_end0':150,'rt_start0':150,'rt_end0':200}])
  j.write_tsv(s/'CANDIDATES.tsv',[{'sample_id':'fixture','match_group':'test','leakage_unit':'test','source_study':'test','family_cluster':'test'}])
  j.write_tsv(s/'SEALED_FILE_SHA256.tsv',j.hashes(r,list(s.rglob('*'))))
  j.write_json(r/'PREFREEZE.json',{'frozen_manifest_sha256':j.sha(f/'FROZEN_FILE_SHA256.tsv'),'truth_seal_manifest_sha256':j.sha(s/'SEALED_FILE_SHA256.tsv')})
  statuses=[]
  for method in j.METHODS:
   dest=o/(method+'.tsv');status='FAILED' if failed and method==j.METHODS[0] else 'SUCCESS'
   if status=='SUCCESS':j.write_tsv(dest,[{'rank':1,'start0':100,'end0':200}])
   statuses.append({'opaque_id':'J2_00001','method':method,'status':status,'output_path':str(dest.relative_to(r))})
  j.write_tsv(o/'METHOD_STATUS.tsv',statuses)
  j.write_tsv(o/'SCORE_FILE_SHA256.tsv',j.hashes(r,list(o.rglob('*'))))
  j.write_json(r/'SCORE_FREEZE.json',{'score_manifest_sha256':j.sha(o/'SCORE_FILE_SHA256.tsv'),'failed_n':int(failed)})
  return r
 def evaluate(self,r):
  with contextlib.redirect_stdout(io.StringIO()):return j.evaluate(types.SimpleNamespace(run=str(r)))
 def test_exact_interval_boundaries(self):
  self.assertFalse(j.overlap(0,100,100,200));self.assertTrue(j.overlap(99,101,100,200))
 def test_complete_failed_denominator_retained(self):
  r=self.make_run(failed=True);self.assertEqual(self.evaluate(r),0)
  metrics=j.read_tsv(r/'evaluation/TARGET_METRICS.tsv');arm=next(x for x in metrics if x['method']==j.METHODS[0])
  self.assertEqual(arm['recovered'],'0');self.assertEqual(arm['status'],'FAILED')
  burden=next(x for x in j.read_tsv(r/'evaluation/ALL_SAMPLE_CANDIDATE_BURDEN.tsv') if x['method']==j.METHODS[0])
  self.assertEqual(burden['candidates_per_mb'],'')
 def test_no_bootstrap_pseudoreplication_for_one_unit(self):
  r=self.make_run();self.evaluate(r)
  self.assertEqual(j.js(r/'evaluation/EVALUATION_STATUS.json')['bootstrap_status'],'DESCRIPTIVE_ONLY_TOO_FEW_INDEPENDENT_POSITIVE_UNITS')
 def test_missing_score_freeze_blocks_truth(self):
  r=self.make_run();(r/'SCORE_FREEZE.json').unlink()
  with self.assertRaises(FileNotFoundError):self.evaluate(r)
  self.assertFalse((r/'evaluation').exists())
 def test_score_mutation_detected(self):
  r=self.make_run();(r/'scores'/f'{j.METHODS[0]}.tsv').write_text('changed\n')
  with self.assertRaises(j.Hold):self.evaluate(r)
 def test_truth_mutation_detected(self):
  r=self.make_run();(r/'sealed/TRUTH.tsv').write_text('changed\n')
  with self.assertRaises(j.Hold):self.evaluate(r)
 def test_incomplete_method_attempts_block_unblind(self):
  r=self.make_run();p=r/'scores/METHOD_STATUS.tsv';rows=j.read_tsv(p);j.write_tsv(p,rows[:-1])
  h=r/'scores/SCORE_FILE_SHA256.tsv';j.write_tsv(h,j.hashes(r,[x for x in (r/'scores').rglob('*') if x!=h]))
  g=j.js(r/'SCORE_FREEZE.json');g['score_manifest_sha256']=j.sha(h);j.write_json(r/'SCORE_FREEZE.json',g)
  with self.assertRaises(j.Hold):self.evaluate(r)
 def test_multicontig_rejected(self):
  p=self.root/'bad.fna';p.write_text('>a\nACGT\n>b\nACGT\n')
  with self.assertRaises(j.Hold):j.fasta(p)
 def test_ambiguous_bases_reported_not_silently_removed(self):
  p=self.root/'N.fna';p.write_text('>a\nACGNT\n');z=j.sequence_summary(p)
  self.assertEqual(z['length_bp'],5);self.assertEqual(z['non_acgt_count'],1)
 def test_empty_external_manifest_never_passes(self):
  cp=self.root/'c.tsv';tp=self.root/'t.tsv';xp=self.root/'x.tsv'
  j.write_tsv(cp,[],j.CANDIDATE_FIELDS);j.write_tsv(tp,[],j.TRUTH_FIELDS);j.write_tsv(xp,[],['key_type','value'])
  cc,tt,obs,errors,notes=j.validate_manifest(cp,tp,xp)
  self.assertTrue(errors)
 def test_secondary_many_loci_is_one_genome(self):
  r=self.make_run();t=r/'sealed/TRUTH.tsv';rows=j.read_tsv(t)
  row=dict(rows[0],label_class='other_structured_rna');row2=dict(row,target_id='fixture_unit_2',truth_start0=300,truth_end0=350)
  j.write_tsv(t,[row,row2]);seal=r/'sealed/SEALED_FILE_SHA256.tsv';j.write_tsv(seal,j.hashes(r,[x for x in (r/'sealed').rglob('*') if x!=seal]))
  pg=j.js(r/'PREFREEZE.json');pg['truth_seal_manifest_sha256']=j.sha(seal);j.write_json(r/'PREFREEZE.json',pg)
  self.evaluate(r);m=j.read_tsv(r/'evaluation/GENOME_AGGREGATE_SUMMARY.tsv')
  self.assertTrue(all(x['native_genome_denominator']=='1' for x in m))
  self.assertEqual(j.js(r/'evaluation/EVALUATION_STATUS.json')['independent_positive_leakage_units'],1)
 def test_background_is_burden_not_negative_accuracy(self):
  r=self.make_run();p=r/'sealed/TRUTH.tsv';rows=j.read_tsv(p);rows[0]['label_class']='background_unknown';j.write_tsv(p,rows)
  seal=r/'sealed/SEALED_FILE_SHA256.tsv';j.write_tsv(seal,j.hashes(r,[x for x in (r/'sealed').rglob('*') if x!=seal]))
  pg=j.js(r/'PREFREEZE.json');pg['truth_seal_manifest_sha256']=j.sha(seal);j.write_json(r/'PREFREEZE.json',pg)
  self.evaluate(r);self.assertEqual(j.read_tsv(r/'evaluation/TARGET_METRICS.tsv'),[])
  burden=j.read_tsv(r/'evaluation/ALL_SAMPLE_CANDIDATE_BURDEN.tsv')
  self.assertEqual(len(burden),len(j.METHODS));self.assertTrue(all(x['label_class']=='background_unknown' for x in burden))
  self.assertTrue(all('precision' not in x and 'false_positive_rate' not in x for x in burden))
 def test_additional_art_target_requires_rt_and_correct_scope(self):
  import random
  rng=random.Random(17);evidence=self.root/'author_fixture.txt';evidence.write_text('synthetic guard-test evidence; not scientific truth')
  cc=[];tt=[]
  for sid,label in [('positive_fixture','positive'),('neighbor_fixture','hard_neighbor')]:
   fa=self.root/(sid+'.fna');fa.write_text('>'+sid+'\n'+''.join(rng.choices('ACGT',k=4096))+'\n');z=j.sequence_summary(fa)
   c={key:'TEST_ONLY' for key in j.CANDIDATE_FIELDS}
   c.update(sample_id=sid,accession=sid,native_parent_accession=sid,genome_path=str(fa),
    fasta_sha256=z['fasta_sha256'],sequence_sha256=z['sequence_sha256'],exposure_status='UNEXPOSED',
    eligibility_status='CONFIRMED',leakage_review_status='REVIEWED_SOURCE_AND_SEQUENCE_GROUPS',
    leakage_review_path=str(evidence),leakage_review_sha256=j.sha(evidence))
   cc.append(c)
   t={key:'' for key in j.TRUTH_FIELDS};t.update(sample_id=sid,target_id=sid+'_1',label_class=label,
    label_scope='author_defined_system_interval',truth_start0=100,truth_end0=500,
    array_start0=100,array_end0=200,rt_start0=250,rt_end0=500,author_truth_url='https://example.invalid/test-only',
    truth_evidence_path=str(evidence),truth_evidence_sha256=j.sha(evidence),
    truth_basis='INDEPENDENT_AUTHOR_ANNOTATION' if label=='positive' else 'INDEPENDENT_NON_ART_ANNOTATION',
    curator='TEST_ONLY',truth_fixed_utc='2026-01-01T00:00:00Z');tt.append(t)
  bad=dict(tt[0],target_id='positive_fixture_2',rt_start0='',rt_end0='',label_scope='wrong_scope');tt.append(bad)
  cp=self.root/'c.tsv';tp=self.root/'t.tsv';xp=self.root/'x.tsv'
  j.write_tsv(cp,cc);j.write_tsv(tp,tt);j.write_tsv(xp,[],['key_type','value'])
  _,_,_,errors,_=j.validate_manifest(cp,tp,xp)
  reasons=[x['reason'] for x in errors]
  self.assertTrue(any('positive_fixture_2 ART-positive requires independent rt' in x for x in reasons))
  self.assertTrue(any('positive_fixture_2 label_scope mismatch' in x for x in reasons))

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
 suite=unittest.defaultTestLoader.loadTestsFromTestCase(Guards);runner=unittest.TextTestRunner(verbosity=2);r=runner.run(suite)
 report={'status':'SOFTWARE_GUARD_TESTS_PASS' if r.wasSuccessful() else 'SOFTWARE_GUARD_TESTS_FAIL',
  'test_n':r.testsRun,'failure_n':len(r.failures),'error_n':len(r.errors),'fixtures':'synthetic temporary files, not scientific validation',
  'new_external_results':False}
 j.write_json(a.output,report);print(json.dumps(report,indent=2));return 0 if r.wasSuccessful() else 1
if __name__=='__main__':sys.exit(main())
