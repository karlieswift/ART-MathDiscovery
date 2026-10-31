from pathlib import Path
import sys,unittest,tempfile,hashlib,json,csv,os,copy
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import stage05b as x

class Tests(unittest.TestCase):
    def test_01_package_and_exact_engine_sources(self):x.package_check()
    def test_02_919_manifest_accounting(self):
        r=x.tsv_read(x.BASE/'config/FROZEN_NEW_GENOME_MANIFEST.tsv')
        self.assertEqual(len(r),919);self.assertEqual(sum(int(a['length_bp']) for a in r),171397507)
    def test_03_all_window_arithmetic(self):
        r=x.tsv_read(x.BASE/'config/FROZEN_NEW_GENOME_MANIFEST.tsv')
        self.assertEqual(sum(x.window_count(int(a['length_bp'])) for a in r),21062200)
        for a in r:self.assertEqual(x.window_count(int(a['length_bp'])),int(a['oriented_windows']))
    def test_04_window_tail(self):
        self.assertEqual(x.starts(161,128),[0,32,33]);self.assertEqual(x.starts(127,128),[])
    def test_05_genome_specific_peak(self):
        r=x.tsv_read(x.BASE/'config/FROZEN_NEW_GENOME_MANIFEST.tsv');p=x.build_plan(r,2,8)
        self.assertEqual(p['max_single_genome_ACNV_memmap_bytes'],1483927712)
        self.assertEqual(p['upper_bound_simultaneous_ACNV_memmap_bytes'],2967855424)
        self.assertEqual(p['cumulative_ACNV_array_bytes'],460672438400)
    def test_06_policy_only_two_arms(self):
        p=x.pol();self.assertEqual(p['families'],['ACNV_FAMILY','GENERIC_RECURRENCE']);self.assertEqual(p['proposal_q_threshold'],.005)
    def test_07_acnv_dimension_real_sequence(self):
        m=x.engine_import();m.init_acnv(str(x.ENGINE/'frozen_sources/ACNV_SOURCE__legacy_acnv_extractor_frozen.py'))
        z=m.acnv_one('ACGTACGCCATG'*12);self.assertEqual(z.shape,(5468,));self.assertEqual(z.dtype,np.float32);self.assertTrue(np.isfinite(z).all())
    def test_08_strict_overlap_not_touch(self):
        m=x.engine_import();d=pd.DataFrame({'start0':[0,10,11],'end0':[10,12,20]});self.assertEqual(m.merge_overlapping(d),[[0,10],[10,20]])
    def test_09_strand_ties_preserve_plus(self):
        m=x.engine_import();d=m.make_windows('ACGT'*100,'SYN',[128],.25);z=m.collapse_oriented_q(d,np.ones(len(d)),'a')
        self.assertTrue((z.a_best_strand=='+').all())
    def test_10_robust_constant(self):
        m=x.engine_import();z,k=m.robust_score(np.ones((20,4),dtype='float32'));self.assertEqual(k,0);self.assertTrue((z==0).all())
    def test_11_robust_nonfinite_fail(self):
        with self.assertRaises(ValueError):x.engine_import().robust_score([[float('nan')]])
    def test_12_manifest_unique_rc_commitments(self):
        r=x.tsv_read(x.BASE/'config/FROZEN_NEW_GENOME_MANIFEST.tsv');self.assertEqual(len({min(a['sequence_sha256'],a['reverse_complement_sha256']) for a in r}),919)
    def test_13_fasta_no_silent_N_repair(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.fa';p.write_text('>A\nACNGT\n')
            with self.assertRaises(RuntimeError):x.read_fasta(p)
    def test_14_fasta_multirecord_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.fa';p.write_text('>A\nACGT\n>B\nACGT\n')
            with self.assertRaises(RuntimeError):x.read_fasta(p)
    def test_15_fasta_no_lowercase_normalization(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.fa';p.write_text('>A\nacgt\n')
            with self.assertRaises(RuntimeError):x.read_fasta(p)
    def test_16_path_escape_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(RuntimeError):x.checked(Path(d),'../secret')
    def test_17_missing_hash_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(RuntimeError):x.file_hash_check(Path(d)/'missing','0'*64)
    def test_18_wrong_hash_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'data';p.write_text('x')
            with self.assertRaises(RuntimeError):x.file_hash_check(p,'0'*64)
    def test_19_historical_fasta_commitments(self):
        for k,v in x.pol()['historical_file_hashes'].items():x.file_hash_check(x.BASE/k,v)
    def test_20_method_signature_stable(self):self.assertEqual(x.method_signature(),x.method_signature())
    def test_21_dataframe_batch_not_statistical_group(self):
        m=x.engine_import();a=np.array([[0],[0],[0],[0],[10],[10],[10],[10]],dtype='float32')
        self.assertFalse(np.array_equal(m.robust_score(a)[0],np.r_[m.robust_score(a[:4])[0],m.robust_score(a[4:])[0]]))
    def test_22_full_frozen_worker_small_synthetic_and_resume(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);f=d/'SYN.fasta';f.write_text('>SYN\n'+'ACGTGCAATTCG'*40+'\n');out=d/'verified'
            z=x.execute_genome('SYN',f,480,d/'job',out,2,d/'log')
            self.assertTrue(z['gate_pass']);self.assertFalse((d/'job/work/SYN/acnv_5468.npy').exists())
            z2=x.execute_genome('SYN',f,480,d/'job',out,1,d/'log');self.assertEqual(z,z2)
            (out/x.BOARDS[1]).write_text('corruption')
            with self.assertRaises(RuntimeError):x.execute_genome('SYN',f,480,d/'job',out,2,d/'log')
    def test_23_no_replacement_source_ledger(self):
        r=x.tsv_read(x.BASE/'config/SOURCE_ACQUISITION_LEDGER.tsv')
        self.assertEqual(len(r),980);self.assertEqual(sum(a['status']=='FROZEN_ELIGIBLE' for a in r),919)
    def test_24_prepublication_date_is_not_new_system(self):
        r=x.tsv_read(x.BASE/'config/FROZEN_NEW_GENOME_MANIFEST.tsv')
        self.assertEqual(sum(a['release_date']>='2026-09-23' for a in r),0)
    def test_25_parent_core_output_is_not_ART_identity(self):
        p=x.pol();self.assertFalse(p['retuning_authorized']);self.assertTrue(p['zero_region_genomes_retained']);self.assertTrue(p['no_rank_topk_survival_cutoff'])
    def test_26_invalid_workers(self):
        with self.assertRaises(RuntimeError):x.build_plan(x.tsv_read(x.BASE/'config/FROZEN_NEW_GENOME_MANIFEST.tsv'),0,8)

if __name__=='__main__':unittest.main()
