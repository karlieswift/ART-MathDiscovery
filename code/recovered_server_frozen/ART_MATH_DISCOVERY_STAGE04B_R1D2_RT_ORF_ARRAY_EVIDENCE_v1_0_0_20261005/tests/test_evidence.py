import unittest,tempfile,json,csv,sys,copy
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import *
from repeat_evidence import *
from biology_tools import parse_prodigal_header,parse_domtbl
from prepare_refs import hmm_meta,split_hmms,unpack_payload
from run_evidence import preflight

class CoordinateTests(unittest.TestCase):
    def test_halfopen_touch_is_not_overlap(self):self.assertFalse(overlaps((1,5),(5,8)))
    def test_touch_distance_zero(self):self.assertEqual(distance((1,5),(5,8)),0)
    def test_fragment_mask_not_outer_span(self):self.assertEqual(covered_bp((12,18),merge_intervals([(0,10),(20,30)])),0)
    def test_mask_union_no_double_count(self):self.assertEqual(covered_bp((0,25),merge_intervals([(2,12),(6,18)])),16)
    def test_reverse_complement(self):self.assertEqual(rc('ACGTTT'),'AAACGT')
    def test_prodigal_negative_header(self):
        x=parse_prodigal_header('NC_123_1 # 7 # 99 # -1 # ID=1_1;partial=00;start_type=GTG')
        self.assertEqual(x[1:4],(6,99,'-'));self.assertEqual(x[4]['start_type'],'GTG')
    def test_prodigal_bad_orientation(self):
        with self.assertRaises(ValueError):parse_prodigal_header('x # 1 # 90 # 0 # ID=1')
    def test_fasta_duplicate_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.fa';p.write_text('>a\nACGT\n>a\nTTTT\n')
            with self.assertRaises(ValueError):fasta_dict(p)
    def test_empty_fasta_record_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.fa';p.write_text('>a\n')
            with self.assertRaises(ValueError):fasta_dict(p)
    def test_manifest_detects_mutation(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);(d/'a').write_text('ok');m=d/'m.tsv';output_manifest(d,[d/'a'],m);assert_manifest(d,m);(d/'a').write_text('no')
            with self.assertRaises(RuntimeError):assert_manifest(d,m)

class RepeatTests(unittest.TestCase):
    def setUp(self):self.p=read_json(ROOT/'config/METHOD_POLICY.json')['repeat_survey']
    def test_three_exact_copies_found(self):
        core='ACGTCAGTGCAT';s='TT'+core+'A'*27+core+'C'*35+core+'GG'
        rr=find_exact_core_groups(s,self.p)
        self.assertTrue(any(r['core']==core and len(r['starts'])==3 for r in rr))
    def test_two_copies_not_three(self):
        core='ACGTCAGTGCAT';rr=find_exact_core_groups(core+'C'*35+core,self.p)
        self.assertFalse(any(r['core']==core for r in rr))
    def test_no_cross_neighborhood_boundary(self):
        core='ACGTCAGTGCAT';s=core+'G'*23+core+'T'*23+core
        self.assertFalse(any(r['core']==core for r in find_exact_core_groups(s[:-5],self.p)))
    def test_max_distance_is_fixed(self):
        core='ACGTCAGTGCAT';s=core+'G'*1100+core+'T'*1100+core
        self.assertFalse(any(r['core']==core for r in find_exact_core_groups(s,self.p)))
    def test_low_complexity_retained_not_promoted(self):
        rr=find_exact_core_groups('A'*250,self.p);self.assertTrue(rr)
        self.assertTrue(all(entropy(r['core'])<1.5 for r in rr))
    def test_entropy(self):self.assertAlmostEqual(entropy('ACGT'*10),2.0)
    def test_nonoverlap(self):
        for r in find_exact_core_groups('ACGT'*100,self.p):
            self.assertTrue(all(b-a>=r['repeat_len'] for a,b in zip(r['starts'],r['starts'][1:])))
    def test_exact_base_identity(self):
        s='TCGATGACCTGA'+'TTTG'*6+'TCGATGACCTGA'+'CACG'*6+'TCGATGACCTGA'
        for r in find_exact_core_groups(s,self.p):self.assertTrue(all(s[p:p+r['repeat_len']]==r['core'] for p in r['starts']))
    def test_run_determinism(self):
        s=('ACGTCAGTGCAT'+'ATCCGG'*7)*5
        self.assertEqual(find_exact_core_groups(s,self.p),find_exact_core_groups(s,self.p))
    def test_invalid_bases_not_a_seed(self):
        self.assertFalse(any('N' in r['core'] for r in find_exact_core_groups('N'*100,self.p)))

class HmmerParserTests(unittest.TestCase):
    def setUp(self):
        self.catalog=[{'name':'RVT_1','accession_version':'PF00078.99','length':100,'role':'RT_CORE'}]
        self.targets={'AA_TEST':'M'+'A'*299}
        self.dom='AA_TEST - 300 RVT_1 PF00078.99 100 1e-20 65.0 0.0 1 1 1e-22 1e-20 64.0 0.0 2 95 21 114 20 115 0.95 test\n'
    def parse(self,s,mode='GA'):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'d';p.write_text('# comment\n'+s);return parse_domtbl(p,self.catalog,self.targets,mode)
    def test_parse_ga(self):
        r=self.parse(self.dom)[0];self.assertEqual(r['evidence_tier'],'RT_CORE_GA');self.assertAlmostEqual(r['hmm_coverage'],.94)
    def test_parse_sensitivity(self):self.assertEqual(self.parse(self.dom,'SENSITIVITY')[0]['evidence_tier'],'RT_CORE_SENSITIVITY')
    def test_weak_not_promoted(self):self.assertEqual(self.parse(self.dom.replace('1e-20','0.5'),'SENSITIVITY')[0]['evidence_tier'],'REPORTED_ONLY')
    def test_aux_not_core(self):
        self.catalog[0]['role']='RT_LIKE_AUXILIARY';self.assertEqual(self.parse(self.dom)[0]['evidence_tier'],'RT_LIKE_AUXILIARY_GA')
    def test_unknown_profile_fails(self):
        with self.assertRaises(ValueError):self.parse(self.dom.replace('RVT_1','wrong'))
    def test_out_of_range_fails(self):
        with self.assertRaises(ValueError):self.parse(self.dom.replace('21 114','21 301'))
    def test_empty_domain_is_valid_no_hit(self):self.assertEqual(self.parse(''),[])
    def test_truncated_row_rejected(self):
        with self.assertRaises(ValueError):self.parse('AA_TEST - 300 RVT_1')
    def test_hmm_missing_GA_rejected(self):
        with self.assertRaises(ValueError):hmm_meta('HMMER3/f\nNAME RVT_1\nACC PF00078.99\nLENG 100\nALPH amino\n//\n')
    def test_hmm_html_rejected(self):
        with self.assertRaises(ValueError):unpack_payload(b'<html>blocked</html>')
    def test_hmm_truncated_rejected(self):
        with self.assertRaises(ValueError):split_hmms('HMMER3/f\nNAME RVT_1\n')
    def test_hmm_header_metadata(self):
        r=hmm_meta('HMMER3/f\nNAME RVT_1\nACC PF00078.99\nLENG 100\nALPH amino\nGA 25.0 25.0;\n//\n');self.assertEqual(r['ga_domain_bits'],25.0)

class RealFrozenInputTests(unittest.TestCase):
    def test_real_full_input_preflight(self):
        j=preflight();self.assertEqual(j['n_candidates'],383);self.assertEqual(j['n_unique_official_AA'],9330)
    def test_genome_only_not_missing_from_candidate_universe(self):
        c=read_tsv(ROOT/'inputs/candidates.tsv');x=[r for r in c if r['sequence_id']=='NC_138483.1'];self.assertEqual(len(x),1);self.assertEqual(x[0]['region_id'],'R0001')
    def test_original_RT_recurrence_is_coding(self):
        reg=read_tsv(ROOT/'inputs/official_protein_registry.tsv');r=next(r for r in reg if r['official_protein_accession']=='YP_010110885.1');self.assertEqual(covered_bp((64960,65088),merge_intervals(parts_from_registry(r))),128)
    def test_full_feature_ledger_count(self):self.assertEqual(len(read_tsv(ROOT/'inputs/all_official_features.tsv')),20804)

if __name__=='__main__':unittest.main(verbosity=2)
