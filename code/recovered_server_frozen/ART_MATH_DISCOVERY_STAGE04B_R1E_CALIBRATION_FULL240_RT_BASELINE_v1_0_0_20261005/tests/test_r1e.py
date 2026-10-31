from pathlib import Path
import unittest, tempfile, sys, shutil, json
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import run_stage as s
from common import *
from repeat_evidence import find_exact_core_groups

class Contracts(unittest.TestCase):
 def test_frozen_denominators(self):
  _,u,c,l=s.package_data();self.assertEqual((len(u),len(c),sum(l.values())),(240,383,16574805))
 def test_candidate_subset(self):
  _,u,c,l=s.package_data();self.assertEqual(len(set(l)-{z['sequence_id'] for z in c}),66)
 def test_controls_are_not_RT_negatives(self):
  r=read_tsv(ROOT/'inputs/controls/RT_REFERENCE_PANEL.tsv');self.assertEqual(len(r),6)
  self.assertEqual(Counter(z['reference_class'] for z in r),{'REGISTERED_ART_ASSOCIATED_RT':3,'REGISTERED_NONARRAY_RT':3})
  self.assertTrue(all(z['RT_negative_control']=='False' for z in r))
 def test_control_aa_and_nt_hashes(self):
  aa=fasta_dict(ROOT/'inputs/controls/RT_REFERENCE_PANEL.faa')
  for r in read_tsv(ROOT/'inputs/controls/RT_REFERENCE_PANEL.tsv'):
   k=r['reference_id'];self.assertEqual(sha_text(aa[k][1]),r['aa_sha256']);self.assertEqual(len(aa[k][1]),int(r['aa_length']))
   self.assertEqual(sha_text(next(fasta(ROOT/'inputs/controls'/f'{k}.cds.fna'))[2]),r['nt_sha256'])
 def test_drt2_registered_sequence_hash(self):
  r=next(r for r in read_tsv(ROOT/'inputs/controls/RT_REFERENCE_PANEL.tsv') if r['name']=='DRT2')
  self.assertEqual(r['aa_sha256'],'13c58226c9ddaeb54ef404c57893843e80f6ae7248f5b5daeef24597d08c0815')
 def test_original_profiles_unchanged(self):
  for m in read_tsv(ROOT/'refs/MODEL_CATALOG.tsv'):self.assertEqual(sha_file(ROOT/m['path']),m['sha256'])
 def test_all174_prior_raw_lock(self):
  rows=read_tsv(ROOT/'config/PRIOR_PRODIGAL_174_RAW_LOCK.tsv');self.assertEqual(len(rows),522)
  self.assertEqual(len({Path(r['path']).parts[2] for r in rows}),174)
 def test_merged_footprints_not_region_counts(self):
  cand=[{'sequence_id':'a','start0':'5','end0':'10'},{'sequence_id':'a','start0':'8','end0':'12'}]
  self.assertEqual(s.masks(cand,{'a':20,'b':20},0),{'a':[[5,12]],'b':[]})
  self.assertEqual(s.masks(cand,{'a':20},10000),{'a':[[0,20]]})
 def test_halfopen_not_boundary_overlap(self):
  self.assertFalse(overlaps((0,5),(5,10)));self.assertTrue(overlaps((0,6),(5,10)))
 def test_three_copy_rule(self):
  pol=read_json(ROOT/'config/METHOD_POLICY.json')['repeat_survey'];motif='ACGTCAGTACGACTGT'
  self.assertFalse(any(s.entropy(r['core'])>=1.5 for r in find_exact_core_groups(motif+'G'*80+motif,pol)))
  self.assertTrue(any(len(r['starts'])>=3 for r in find_exact_core_groups(motif+'G'*80+motif+'C'*80+motif,pol)))
 def test_position_null_reproducible(self):
  cand=[{'sequence_id':'a','start0':'10000','end0':'11000'}];ll={'a':100000,'b':50000}
  rt=[{'sequence_id':'a','start0':10100,'end0':10200},{'sequence_id':'b','start0':100,'end0':200}]
  r,j=s.position_null(cand,ll,rt,20,123);r2,j2=s.position_null(cand,ll,rt,20,123)
  self.assertEqual((r,j),(r2,j2));self.assertEqual(j['observed_captured_predicted_core_loci'],1)
  self.assertTrue(all(z['captured_predicted_core_loci']<=1 for z in r))
 def test_position_null_zero_loci(self):
  r,j=s.position_null([{'sequence_id':'a','start0':0,'end0':10}],{'a':100},[],10)
  self.assertIsNone(j['enrichment_vs_random_mean']);self.assertTrue(j['no_informative_RT_reference_loci'])
 def test_lock_refuses_changes(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'LOCK.json';s.lock_json(p,{'x':1});s.lock_json(p,{'x':1})
   with self.assertRaises(RuntimeError):s.lock_json(p,{'x':2})
 def test_prediction_minus_orientation(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);p=root/'work/prodigal/ACC';p.mkdir(parents=True);(root/'inputs/genomes').mkdir(parents=True)
   write_fasta(root/'inputs/genomes/ACC.fna',[('ACC','TTATTTCAT')])
   hdr='ACC_1 # 1 # 9 # -1 # partial=00;start_type=ATG'
   (p/'proteins.faa').write_text('>'+hdr+'\nMK*\n');(p/'cds.fna').write_text('>'+hdr+'\nATGAAATAA\n')
   (p/'genes.gff').write_text('# transl_table=11\nACC\tProdigal\tCDS\t1\t9\t.\t-\t0\tID=1_1\n')
   rows,aa=s.parse_prediction(root,'ACC');self.assertEqual(aa['ACC__PRODIGAL_00001'],'MK');self.assertEqual(rows[0]['strand'],'-')
   (p/'cds.fna').write_text('>'+hdr+'\nATGAAATAG\n')
   with self.assertRaises(RuntimeError):s.parse_prediction(root,'ACC')
 def test_unsafe_source_path_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);a04=root/'a04';(a04/'outputs').mkdir(parents=True);(a04/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv').write_text('x')
   pol={'stage04a_manifest_sha256':sha_file(a04/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv')}
   with patch.object(s,'package_data',return_value=(pol,[{'frozen_fasta':'../bad','sequence_id':'ACC'}],[],{'ACC':9})):
    with self.assertRaises(RuntimeError):s.preflight_genomes(a04,root)
 def test_domain_parser_checks_model_both_ga_scores(self):
  cat=[{'name':'RVT_1','length':'200','accession_version':'PF00078.33','role':'RT_CORE','ga_sequence_bits':'29.6','ga_domain_bits':'29.6'}]
  # target, targetacc,len, query, queryacc,M, fullE,fullS,bias,num,of,cE,iE,domS,bias,hmmfrom,to,alifrom,to,envfrom,to,acc
  v=['AA_a','-','346','RVT_1','PF00078.33','200','1e-14','53.4','0','1','1','1e-14','2e-14','52.9','0','1','186','80','270','78','272','0.97','desc']
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'dom';p.write_text(' '.join(v)+'\n');r=s.parse_domains(p,cat,{'AA_a':'A'*346},'GA')
   self.assertEqual(r[0]['tier'],'CORE_GA');self.assertEqual(r[0]['hmm_coverage'],.93)
   v[13]='21.0';p.write_text(' '.join(v)+'\n')
   with self.assertRaises(RuntimeError):s.parse_domains(p,cat,{'AA_a':'A'*346},'GA')
 def test_auxiliary_not_core(self):
  cat=[{'name':'RVT_3','length':'123','accession_version':'PF13456.13','role':'RT_LIKE_AUXILIARY','ga_sequence_bits':'21.6','ga_domain_bits':'21.6'}]
  v=['AA_a','-','188','RVT_3','PF13456.13','123','1e-7','30','0','1','1','1e-7','1e-6','30','0','1','73','1','80','1','90','0.97','d']
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'dom';p.write_text(' '.join(v)+'\n');r=s.parse_domains(p,cat,{'AA_a':'A'*188},'GA');self.assertEqual(r[0]['tier'],'AUX_GA')
 def test_hmmer_command_fixed_universe_and_resume(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);(root/'outputs').mkdir();(root/'refs').mkdir();(root/'logs').mkdir()
   shutil.copyfile(ROOT/'refs/RT_PANEL.hmm',root/'refs/RT_PANEL.hmm');shutil.copyfile(ROOT/'refs/MODEL_CATALOG.tsv',root/'refs/MODEL_CATALOG.tsv')
   write_fasta(root/'outputs/SCAN_TARGETS.faa',[('AA_fake','AAAA')]);calls=[]
   def fake_run(cmd,log):
    calls.append([str(x) for x in cmd]);self.assertEqual(cmd[cmd.index('-Z')+1],'100000')
    rep=Path(cmd[cmd.index('-o')+1]);rep.write_text(''.join('Query: '+x+' [M=200]\nTarget sequences: 1\n' for x in ['RVT_1','RVT_2','RVT_3'])+'[ok]\n')
    Path(cmd[cmd.index('--domtblout')+1]).write_text('# no hits\n');Path(cmd[cmd.index('--tblout')+1]).write_text('# no hits\n')
   with patch.object(s,'run_cmd',side_effect=fake_run):
    s.hmmer(root,{'AA_fake':'AAAA'},1,'testlock');s.hmmer(root,{'AA_fake':'AAAA'},1,'testlock')
   self.assertEqual(len(calls),2);self.assertIn('--cut_ga',calls[0]);self.assertIn('--domE',calls[1])
 def test_actual_source_repeat_precompute_reproducible(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);shutil.copytree(ROOT/'inputs/controls',root/'inputs/controls');(root/'config').mkdir();(root/'outputs').mkdir()
   shutil.copyfile(ROOT/'config/METHOD_POLICY.json',root/'config/METHOD_POLICY.json')
   a=s.repeat_calibration(root);h=sha_file(root/'outputs/REFERENCE_REPEAT_GROUPS.tsv');b=s.repeat_calibration(root)
   self.assertEqual(a,b);self.assertEqual(sha_file(root/'outputs/REFERENCE_REPEAT_GROUPS.tsv'),h)
from collections import Counter
if __name__=='__main__':unittest.main()
