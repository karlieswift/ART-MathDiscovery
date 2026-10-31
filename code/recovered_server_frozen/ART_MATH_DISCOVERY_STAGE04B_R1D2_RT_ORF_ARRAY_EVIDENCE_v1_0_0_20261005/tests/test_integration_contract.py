"""Real input cardinality / synthetic domain join regression, no real scan claims."""
import unittest,tempfile,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import *
from biology_tools import build_targets
from integrate_evidence import integrate

class IntegrationTests(unittest.TestCase):
    def make_root(self,td):
        t=Path(td);(t/'inputs').symlink_to(ROOT/'inputs',target_is_directory=True);(t/'outputs').mkdir();return t
    def test_full_official_dedup_and_record_fanout(self):
        with tempfile.TemporaryDirectory() as td:
            t=self.make_root(td);targets,fan=build_targets(t,[],{})
            self.assertEqual((len(targets),len(fan)),(9330,9576))
            self.assertEqual(sum(x['feature_type']=='MATURE_PEPTIDE' for x in fan),2)
            self.assertEqual(sum(len(x['candidate_keys'].split(';')) for x in fan),11708)
    def test_target_freeze_resume_is_immutable(self):
        with tempfile.TemporaryDirectory() as td:
            t=self.make_root(td);build_targets(t,[],{});lp=t/'outputs/TARGET_UNIVERSE_LOCK.json';b=lp.read_bytes();mtime=lp.stat().st_mtime_ns
            build_targets(t,[],{});self.assertEqual(lp.read_bytes(),b);self.assertEqual(lp.stat().st_mtime_ns,mtime)
    def test_target_tamper_is_not_silently_repaired(self):
        with tempfile.TemporaryDirectory() as td:
            t=self.make_root(td);build_targets(t,[],{});p=t/'outputs/RT_SCAN_TARGETS.faa';p.write_text('changed\n')
            with self.assertRaises(RuntimeError):build_targets(t,[],{})
            self.assertEqual(p.read_text(),'changed\n')
    def test_all383_original_fields_survive_synthetic_domain_join(self):
        with tempfile.TemporaryDirectory() as td:
            t=self.make_root(td);targets,fan=build_targets(t,[],{});c=read_tsv(t/'inputs/candidates.tsv')
            windows=[]
            for r in c:
                windows.append({'candidate_key':r['sequence_id']+'::'+r['region_id'],'official_annotation_available':int(r['sequence_id']!='NC_138483.1'),'n_exact_core_repeat_groups':0,'n_high_complexity_repeat_groups':0,'n_high_complexity_full_span_outside_official_cds':0,'n_high_complexity_full_span_outside_official_and_predicted_cds':0,'recurrence_official_cds_overlap_bp':0,'recurrence_predicted_cds_overlap_bp':0})
            tid=fan[0]['target_id']
            g=integrate(t,fan,[{'target_id':tid,'evidence_tier':'RT_LIKE_AUXILIARY_SENSITIVITY'}],[],windows,[])
            out=read_tsv(t/'outputs/ALL_383_RT_ORF_REPEAT_EVIDENCE.tsv')
            self.assertEqual(len(out),383)
            for old,new in zip(c,out):self.assertTrue(all(new[k]==v for k,v in old.items()))
            labels={r['evidence_label'] for r in out if r['candidate_key'] in fan[0]['candidate_keys'].split(';')}
            self.assertEqual(labels,{'RT_LIKE_AUXILIARY_SENSITIVITY_ONLY'})
            self.assertFalse(g['ART_identity_classification_complete'])

if __name__=='__main__':unittest.main()
