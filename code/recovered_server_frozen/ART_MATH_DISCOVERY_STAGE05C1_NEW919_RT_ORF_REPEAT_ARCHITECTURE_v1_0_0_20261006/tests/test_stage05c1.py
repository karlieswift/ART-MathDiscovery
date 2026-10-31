import copy,json,math,random,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import stage05c1 as s
from approximate_repeats import find_approximate_groups,hamming_at_most,greedy_components
from repeat_evidence import find_exact_core_groups,entropy

class Stage05C1Tests(unittest.TestCase):
    def setUp(self):self.pol=s.policy()
    def dna(self,n,seed=11):
        r=random.Random(seed);return ''.join(r.choice('ACGT') for _ in range(n))
    def cand(self):
        return dict(sequence_id='TEST.1',region_id='AR0001',start0='10',end0='140',recurrence_best_window_start0='0',recurrence_best_window_end0='200',recurrence_best_k='8',recurrence_best_lag_bp='40')
    def test_01_frozen_candidate_membership(self):
        c=s.read_tsv(s.ROOT/'config/STAGE05B_DUAL_ARM_CANDIDATES.tsv');self.assertEqual(len(c),3674);self.assertEqual(len({r['sequence_id'] for r in c}),828)
    def test_02_frozen_manifest_hash(self):s.verify_file(s.ROOT/'config/FROZEN_NEW_GENOME_MANIFEST.tsv',self.pol['source_manifest_sha256'])
    def test_03_candidate_hash(self):s.verify_file(s.ROOT/'config/STAGE05B_DUAL_ARM_CANDIDATES.tsv',self.pol['source_candidate_sha256'])
    def test_04_profile_bytes(self):
        for r in s.read_tsv(s.ROOT/'refs/MODEL_CATALOG.tsv'):s.verify_file(s.ROOT/'refs'/(r['accession']+'.hmm'),r['sha256'])
    def test_05_reference_protein_count(self):self.assertEqual(len(s.fasta_dict(s.ROOT/'inputs/controls/RT_REFERENCE_PANEL.faa')),6)
    def test_06_actual_reference_no_negative_label(self):
        self.assertTrue(all(r['RT_negative_control']=='False' for r in s.read_tsv(s.ROOT/'inputs/controls/RT_REFERENCE_PANEL.tsv')))
    def test_07_rc(self):self.assertEqual(s.rc('ATGCCA'),'TGGCAT')
    def test_08_half_open_touch(self):self.assertFalse(s.overlaps((1,5),(5,9)));self.assertEqual(s.distance((1,5),(5,9)),0)
    def test_09_union_coverage(self):self.assertEqual(s.covered_bp((2,15),s.merge_intervals([(1,7),(5,10),(12,20)])),11)
    def test_10_path_escape(self):
        with self.assertRaises(RuntimeError):s.safe_path(s.ROOT,'../escape')
    def test_11_parse_plus_header(self):
        _,a,b,st,_=s.parse_pred_header('X_1 # 11 # 52 # 1 # ID=1_1;partial=00;');self.assertEqual((a,b,st),(10,52,'+'))
    def test_12_parse_minus_header(self):
        _,a,b,st,at=s.parse_pred_header('X_2 # 21 # 120 # -1 # ID=1_2;partial=11;');self.assertEqual((a,b,st,at['partial']),(20,120,'-','11'))
    def test_13_bad_header_rejected(self):
        with self.assertRaises(RuntimeError):s.parse_pred_header('X # 1 # 9 # 0 # ID=x')
    def test_14_hamming_boundary(self):self.assertEqual(hamming_at_most('ACGT','ACGA',1),1);self.assertIsNone(hamming_at_most('ACGT','AAAA',1))
    def test_15_hamming_length(self):self.assertIsNone(hamming_at_most('AC','ACG',3))
    def test_16_nonoverlap_components(self):self.assertEqual(greedy_components([0,1,20,80],12,40),[[0,20],[80]])
    def test_17_exact_three_copies(self):
        m='ACGTTGCATGCAACGT';z=m+self.dna(80,1)+m+self.dna(80,2)+m
        self.assertTrue(any(r['repeat_len']>=16 and len(r['starts'])>=3 for r in find_exact_core_groups(z,self.pol['repeat_survey'])))
    def test_18_exact_two_not_three(self):
        m='ACGTTGCATGCAACGT';z=m+self.dna(80)+m;self.assertEqual(find_exact_core_groups(z,self.pol['repeat_survey']),[])
    def test_19_entropy_low_not_high(self):self.assertEqual(entropy('A'*32),0);self.assertEqual(entropy('ACGT'*8),2)
    def test_20_approx_substitutions(self):
        m='ACGTTGCATGCAACGTACGTCAGT';b=list(m);b[10]='A' if b[10]!='A' else 'C';c=list(m);c[19]='A' if c[19]!='A' else 'C';z=m+self.dna(80,3)+''.join(b)+self.dna(80,4)+''.join(c)
        r=find_approximate_groups(z,self.pol['approximate_repeat_sensitivity']);self.assertTrue(any(v['repeat_len']==24 and len(v['starts'])>=3 and max(v['mismatches'])>0 for v in r))
    def test_21_approx_emitted_boundaries(self):
        m='ACGTTGCATGCAACGT';z=m+self.dna(80,5)+m+self.dna(80,6)+m
        for r in find_approximate_groups(z,self.pol['approximate_repeat_sensitivity']):
            self.assertTrue(all(hamming_at_most(r['core'],z[i:i+r['repeat_len']],math.floor(r['repeat_len']*.2)) is not None for i in r['starts']))
    def test_22_approx_two_no_array(self):
        z='ACGTTGCATGCAACGT'+self.dna(90,9)+'ACGTTGCATGCAACGT';self.assertFalse(find_approximate_groups(z,self.pol['approximate_repeat_sensitivity']))
    def test_23_short_sequence(self):self.assertEqual(find_approximate_groups('ACGT',self.pol['approximate_repeat_sensitivity']),[])
    def test_24_approx_no_RT_input(self):self.assertNotIn('RT',find_approximate_groups.__code__.co_varnames)
    def test_25_witness_matches_not_copies(self):
        z=list(self.dna(200,13));mot='ATGAAAATCTGGTCCGG';z[10:27]=mot;z[94:111]=mot;c=self.cand();c['recurrence_best_lag_bp']='84';w=s.recurrence_witness(''.join(z),c,[(0,200)])
        blocks=json.loads(w['blocks_json']);self.assertTrue(any(b['block_length']>=17 and b['n_overlapping_kmer_matches']>=10 for b in blocks));self.assertFalse(w['matching_kmer_count_is_repeat_copy_count']);self.assertEqual(w['window_predicted_CDS_overlap_bp'],200)
    def test_26_witness_empty(self):
        z=self.dna(200,77);w=s.recurrence_witness(z,self.cand(),[]);self.assertEqual(w['n_matching_kmer_pairs'],0)
    def test_27_coding_span_vs_copy_overlap(self):
        z='ACGT'*60;g={'repeat_len':16,'starts':(0,40,80),'core':'ACGT'*4,'seed_count':1};r=s.annotate_repeat(z,self.cand(),0,g,'EXACT_CORE_PRIMARY_DESCRIPTION',1,[(20,30)])
        self.assertEqual(r['predicted_CDS_overlap_copies_bp'],0);self.assertEqual(r['predicted_CDS_overlap_span_bp'],10);self.assertFalse(r['full_span_outside_predicted_CDS'])
    def test_28_bad_repeat_sequence(self):
        g={'repeat_len':16,'starts':(0,40,80),'core':'A'*16,'seed_count':1}
        with self.assertRaises(RuntimeError):s.annotate_repeat('ACGT'*60,self.cand(),0,g,'EXACT_CORE_PRIMARY_DESCRIPTION',1,[])
    def test_29_checkpoint_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'x').write_text('A');s.write_json(p/'COMPLETE.json',{'fingerprint':{'x':1},'files':{'x':s.sha_file(p/'x')}});self.assertTrue(s.checkpoint(p,{'x':1},['x']));(p/'x').write_text('B')
            with self.assertRaises(RuntimeError):s.checkpoint(p,{'x':1},['x'])
    def test_30_checkpoint_signature(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);s.write_json(p/'COMPLETE.json',{'fingerprint':{'x':1},'files':{}})
            with self.assertRaises(RuntimeError):s.checkpoint(p,{'x':2},[])
    def domain(self,seqscore=40,domscore=40,ie='1e-8',hf=1,ht=199):
        return f'AA_test - 346 RVT_1 PF00078.33 200 1e-9 {seqscore} 0 1 1 1e-10 {ie} {domscore} 0 {hf} {ht} 80 270 70 280 0.90 -\n'
    def domain_parse(self,line,mode):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';p.write_text(line);return s.parse_domains(p,{'AA_test':'A'*346},s.read_tsv(s.ROOT/'refs/MODEL_CATALOG.tsv'),mode)
    def test_31_GA_both_thresholds(self):
        self.assertEqual(self.domain_parse(self.domain(),'GA')[0]['evidence_tier'],'RT_CORE_GA')
        with self.assertRaises(RuntimeError):self.domain_parse(self.domain(domscore=20),'GA')
    def test_32_sensitivity_coverage(self):
        self.assertEqual(self.domain_parse(self.domain(hf=1,ht=40),'SENSITIVITY')[0]['evidence_tier'],'REPORTED_ONLY')
    def test_33_sensitivity_threshold(self):self.assertEqual(self.domain_parse(self.domain(ie='.01'),'SENSITIVITY')[0]['evidence_tier'],'REPORTED_ONLY')
    def test_34_unknown_profile(self):
        with self.assertRaises(RuntimeError):self.domain_parse(self.domain().replace('PF00078.33','PF00078.34'),'GA')
    def test_35_status_no_doublecount(self):
        r=self.domain_parse(self.domain(),'GA')+self.domain_parse(self.domain(),'SENSITIVITY');flags=s.statuses(r);self.assertEqual(flags['AA_test'],{'RT_CORE_GA','RT_CORE_SENSITIVITY'})
    def test_36_same_prediction_two_candidate_links(self):
        seq='ATGAAATAA';desc='TEST.1_1 # 1 # 9 # 1 # ID=1_1;partial=00;start_type=ATG;'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'proteins.faa').write_text('>'+desc+'\nMK*\n');(p/'cds.fna').write_text('>'+desc+'\n'+seq+'\n');(p/'genes.gff').write_text('# Model Data: transl_table=11;\nTEST.1\tProdigal\tCDS\t1\t9\t.\t+\t0\tID=1_1\n')
            rows,aas=s.parse_prediction(p,'TEST.1',seq,[{'region_id':'AR1','start0':'1','end0':'3'},{'region_id':'AR2','start0':'4','end0':'7'}]);self.assertEqual(len(rows),1);self.assertEqual(rows[0]['candidate_keys'],'TEST.1::AR1;TEST.1::AR2');self.assertEqual(next(iter(aas.values())),'MK')
    def test_37_internal_stop_not_replaced(self):
        seq='ATGAAATAA';desc='TEST.1_1 # 1 # 9 # 1 # ID=1_1;partial=00;'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'proteins.faa').write_text('>'+desc+'\nM*K\n');(p/'cds.fna').write_text('>'+desc+'\n'+seq+'\n');(p/'genes.gff').write_text('# transl_table=11;\nTEST.1\tP\tCDS\t1\t9\t.\t+\t0\tID=1_1\n')
            with self.assertRaises(RuntimeError):s.parse_prediction(p,'TEST.1',seq,[])
    def test_38_pair_direction(self):
        l={'start0':100,'end0':200,'strand':'-','record_key':'P','profile_evidence':'RT_CORE_GA','exact_reference_AA_ids':''};r={'start0':210,'end0':280,'repeat_group_id':'R','assay':'EXACT','n_copies':3,'min_identity_to_template':1,'high_complexity':1,'full_span_outside_predicted_CDS':1,'overlap_candidate_bp':0};p=s.pair_context(self.cand(),l,r);self.assertEqual(p['repeat_relative_to_RT_coding_direction'],'UPSTREAM_RELATIVE_TO_PREDICTED_RT');self.assertEqual(p['distance_RT_to_repeat_bp'],10)
    def test_39_duplicate_fasta_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';p.write_text('>x\nACG\n>x\nCGT\n')
            with self.assertRaises(ValueError):s.fasta_dict(p)
    def test_40_original_primary_untouched_policy(self):self.assertFalse(self.pol['mathematical_rescoring']);self.assertFalse(self.pol['UG27_touched']);self.assertFalse(self.pol['negative_calls_authorized'])
    def test_41_reference_exact_counts(self):
        for r in s.read_tsv(s.ROOT/'inputs/controls/REPEAT_REFERENCE_CONTEXT.tsv'):
            seq=next(s.fasta(s.ROOT/r['path']))[2];g=find_exact_core_groups(seq,self.pol['repeat_survey']);self.assertEqual(len(g),{'REF_Madawaska':2,'REF_Machias':0,'REF_SA1':4}[r['reference_id']])
    def test_42_complete_scope_all919(self):self.assertEqual(self.pol['n_genomes'],919);self.assertTrue(self.pol['no_RT_gate_for_repeat_scan']);self.assertTrue(self.pol['no_repeat_gate_for_RT_locus_preservation'])

    def test_43_context_hash_is_nucleotide_content(self):
        for r in s.read_tsv(s.ROOT/'inputs/controls/REPEAT_REFERENCE_CONTEXT.tsv'):
            seq=next(s.fasta(s.ROOT/r['path']))[2]
            self.assertEqual(s.sha_text(seq),r['context_sha256'])
    def test_44_nohit_candidate_retained(self):
        c=self.cand();key=c['sequence_id']+'::'+c['region_id'];w=s.recurrence_witness(self.dna(200,77),c,[]);ctx={'candidate_key':key,'left_clipped':1,'right_clipped':0}
        with tempfile.TemporaryDirectory() as d:
            with patch.object(s,'ROOT',Path(d)):
                out,pairs=s.integrate([], [c], [], [], [], [w], [ctx])
                self.assertEqual(len(out),1);self.assertEqual(pairs,[]);self.assertEqual(out[0]['evidence_label'],'NO_JOINT_EVIDENCE_IN_THESE_LIMITED_ASSAYS');self.assertFalse(out[0]['negative_call_authorized'])
                self.assertTrue(all(out[0][k]==v for k,v in c.items()))
    def test_45_approx_policy_not_exhaustive(self):
        self.assertFalse(self.pol['approximate_repeat_sensitivity']['indels_supported']);self.assertFalse(self.pol['approximate_repeat_sensitivity']['RNA_structure_computed']);self.assertFalse(self.pol['approximate_repeat_sensitivity']['candidate_exclusion_on_absence'])

    def test_46_spacer_diversity_and_empty_spacers(self):
        z='ACGT'*60;g={'repeat_len':16,'starts':(0,16,40),'core':'ACGT'*4,'seed_count':1};r=s.annotate_repeat(z,self.cand(),0,g,'EXACT_CORE_PRIMARY_DESCRIPTION',1,[])
        self.assertEqual(r['n_spacers'],2);self.assertEqual(r['n_empty_spacers'],1);self.assertEqual(r['n_distinct_spacer_sequences'],2)
        self.assertEqual(json.loads(r['spacer_lengths_json']),[0,8])

if __name__=='__main__':unittest.main(verbosity=2)
