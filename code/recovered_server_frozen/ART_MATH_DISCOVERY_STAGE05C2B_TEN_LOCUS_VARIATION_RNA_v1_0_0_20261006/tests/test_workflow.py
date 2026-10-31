import unittest,sys,tempfile,random,json,csv
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import *
from prepare import align_metrics,ORDER,CLOSE
from fold import shuffled_dinucs,dinucs,parse_output
from author_sources import extract_pdf

class WorkflowTests(unittest.TestCase):
 def test_count_tree(self):
  x=read_json(ROOT/'precomputed/COUNT_DEFINITIONS.json')
  self.assertEqual(x['whole_genome_RT_core_GA_loci'],x['math_neighborhood_captured']+x['math_neighborhood_not_captured'])
  self.assertEqual(7,5+1+1);self.assertEqual(6,5+1)
 def test_ten_loci_unique(self):
  x=read_tsv(ROOT/'precomputed/TEN_LOCUS_UNAMBIGUOUS_LEDGER.tsv');self.assertEqual(len(x),10);self.assertEqual(len({r['sequence_id'] for r in x}),10)
 def test_two_definitions_of_five(self):
  x=read_tsv(ROOT/'precomputed/TEN_LOCUS_UNAMBIGUOUS_LEDGER.tsv')
  exact={r['sequence_id'] for r in x if r['C1_exact_plus_RT_label']=='True'}
  close={r['sequence_id'] for r in x if r['C2A_close_three_component_reference_relative']=='True' and r['math_neighborhood_captured']=='True'}
  self.assertEqual(len(exact),5);self.assertEqual(len(close),5);self.assertEqual(exact-close,{'NC_016163.1'});self.assertEqual(close-exact,{'OR455461.1'})
 def test_uncaptured_is_not_eleventh(self):
  x={r['sequence_id']:r for r in read_tsv(ROOT/'precomputed/TEN_LOCUS_UNAMBIGUOUS_LEDGER.tsv')}
  self.assertEqual(x['OR770614.1']['math_neighborhood_captured'],'False');self.assertEqual(x['OR770614.1']['C2A_close_three_component_reference_relative'],'True');self.assertEqual(len(x),10)
 def test_remaining_3674_not_discarded(self):
  self.assertTrue(read_json(ROOT/'config/POLICY.json')['invariants']['all_3674_candidates_retained_upstream'])
 def test_global_exact(self):
  m,v=align_metrics('ACGTACGT','ACGTACGT',False);self.assertEqual(m['identity_all_columns_pct'],100);self.assertEqual(v,[])
 def test_global_substitution(self):
  m,v=align_metrics('A'*20+'C'+'G'*20,'A'*20+'T'+'G'*20,False);self.assertEqual(m['mismatch_columns'],1)
 def test_global_insertion(self):
  m,v=align_metrics('A'*20+'G'*20,'A'*20+'TTT'+'G'*20,False);self.assertEqual(m['gap_columns'],3)
 def test_pairwise_count(self):
  x=read_tsv(ROOT/'precomputed/SIX_CLOSE_LOCUS_GLOBAL_COMPONENT_COMPARISON.tsv');self.assertEqual(len(x),45);self.assertEqual(Counter(r['component'] for r in x),{'RT_PROTEIN':15,'PARTNER_LIKE_PROTEIN':15,'FULL_UPSTREAM_DNA':15})
 def test_upstream_variants_reconstruct_target_independently(self):
  seq=fasta(ROOT/'inputs/evidence/TEN_RT_ORIENTED_UPSTREAM_INTERGENIC.fna');vs=read_tsv(ROOT/'precomputed/SIX_CLOSE_UPSTREAM_VARIANT_COLUMNS.tsv')
  pairs={(r['a'],r['b']) for r in vs}
  self.assertEqual(len(pairs),15)
  for a,b in pairs:
   changes=[r for r in vs if (r['a'],r['b'])==(a,b)];p=q=0;out=[]
   for r in changes:
    x,u,y,v=(int(r[k]) for k in ['a_start','a_end','b_start','b_end'])
    self.assertEqual(seq[a][p:x],seq[b][q:y]);self.assertEqual(seq[a][x:u],r['a_seq']);self.assertEqual(seq[b][y:v],r['b_seq']);out.append(seq[a][p:x]);out.append(r['b_seq']);p,q=u,v
   self.assertEqual(seq[a][p:],seq[b][q:]);out.append(seq[a][p:]);self.assertEqual(''.join(out),seq[b])
 def test_merge_and_complement(self):
  self.assertEqual(complement(0,20,[(2,7),(4,8),(10,13),(-8,1),(19,40)]),[(1,2),(8,10),(13,19)])
 def test_window_tail(self):
  x=windows(50,309,[120]);self.assertEqual(x[-1][:2],(189,309));self.assertEqual(len(set(x)),len(x))
 def test_short_windows(self):
  self.assertEqual(windows(0,1),[]);self.assertEqual(windows(5,100),[(5,100,'SHORT_INTERGENIC')])
 def test_reverse_complement(self):self.assertEqual(rc('AACCGT'),'ACGGTT')
 def test_RNA_freeze_count_hashes(self):
  p=ROOT/'precomputed';f=read_json(p/'INPUT_FREEZE.json');s=fasta(p/'RNA_TARGETS.fasta');r=read_tsv(p/'RNA_OCCURRENCE_REGISTRY.tsv')
  self.assertEqual(len(s),f['n_unique_RNA_targets']);self.assertEqual(len(r),f['n_RNA_occurrences']);self.assertEqual(sha_file(p/'RNA_TARGETS.fasta'),f['RNA_targets_sha256'])
 def test_occurrence_real_nt_identity_both_directions(self):
  e=ROOT/'inputs/evidence';f=ROOT/'inputs/focused_sequences';meta={r['sequence_id']:r for r in read_tsv(e/'TEN_RT_CONTEXT_SEQUENCE_MANIFEST.tsv')};refmeta={r['reference_id']:r for r in read_tsv(ROOT/'inputs/reference_sources/REPEAT_REFERENCE_CONTEXT.tsv')};refs=fasta(e/'REGISTERED_ARRAY_CORE_REFERENCES.fna');ctx={k:next(iter(fasta(f/r['file']).values())) for k,r in meta.items()};targets=fasta(ROOT/'precomputed/RNA_TARGETS.fasta')
  for r in read_tsv(ROOT/'precomputed/RNA_OCCURRENCE_REGISTRY.tsv'):
   a,b=int(r['start0']),int(r['end0']);k=r['sequence_id']
   if r['registered_reference']=='True':off=int(refmeta[k]['registered_core_start0']);nt=refs[k][a-off:b-off]
   else:off=int(meta[k]['start0']);nt=ctx[k][a-off:b-off]
   if r['genomic_strand_hypothesis']=='-':nt=rc(nt)
   self.assertEqual(nt.replace('T','U'),targets[r['target_id']])
 def test_both_RNA_orientations_retained(self):
  groups={}
  for r in read_tsv(ROOT/'precomputed/RNA_OCCURRENCE_REGISTRY.tsv'):groups.setdefault((r['sequence_id'],r['start0'],r['end0']),set()).add(r['genomic_strand_hypothesis'])
  self.assertTrue(all(s=={'+','-'} for s in groups.values()))
 def test_all_ten_plus_three_have_RNA(self):
  x=read_tsv(ROOT/'precomputed/RNA_INPUT_COUNTS.tsv');self.assertEqual(len(x),13);self.assertTrue(all(int(r['RNA_occurrences'])>0 for r in x))
 def test_dinucleotide_preservation(self):
  g=random.Random(1)
  for length in [1,2,60,120,480,1285]:
   s=''.join(g.choice('ACGU') for _ in range(length))
   for seed in range(10):
    t=shuffled_dinucs(s,seed);self.assertEqual(dinucs(t),dinucs(s));self.assertEqual(Counter(t),Counter(s));self.assertEqual(t[0],s[0]);self.assertEqual(t[-1],s[-1])
 def test_dinucleotide_deterministic_and_degenerate(self):
  self.assertEqual(shuffled_dinucs('ACGAUCGUACG',42),shuffled_dinucs('ACGAUCGUACG',42));self.assertEqual(shuffled_dinucs('A'*120,7),'A'*120)
 def test_RNA_parser(self):
  x=parse_output('>X\nGCAU\n(..) (-0.20)\n',{'X':'GCAU'});self.assertEqual(x['X']['MFE_kcal_mol'],-.2)
 def test_RNA_parser_fail_closed(self):
  for raw in ['>X\nGCAA\n(..) (-1.0)\n','>X\nGCAU\n(... (-1.0)\n','>X\nGCAU\n(..) (-1.0)\n>X\nGCAU\n(..) (-1.0)\n','garbage']:
   with self.assertRaises(ValueError):parse_output(raw,{'X':'GCAU'})
 def test_RNA_no_p_value_or_novelty(self):
  p=read_json(ROOT/'config/POLICY.json');self.assertIn('no p-values',p['RNAfold']['null_interpretation']);self.assertFalse(p['invariants']['full_catalogue_novelty_audit_complete'])
 def test_source_PDF_inventory_not_membership(self):
  import fitz
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);d=fitz.open();page=d.new_page();page.insert_text((30,30),'Test fixture: OR770614.1 and Data availability');page.insert_link({'kind':fitz.LINK_URI,'from':fitz.Rect(20,20,100,50),'uri':'https://example.org/data'});d.save(p/'test.pdf');d.close()
   info=extract_pdf(p/'test.pdf',p/'outputs/author',['OR770614.1','OR455461.1']);self.assertFalse(info['full_family_catalogue_acquired_and_validated']);x=read_tsv(p/'outputs/author/TEN_LOCUS_AUTHOR_DOCUMENT_STATUS.tsv');self.assertEqual(x[0]['exact_version_mentioned'],'True');self.assertEqual(x[1]['novelty_confirmed'],'False')
 def test_false_PDF_rejected(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'x.pdf';p.write_text('<html>blocked</html>')
   with self.assertRaises(ValueError):extract_pdf(p,Path(td)/'out',ORDER)
 def test_manifest_corruption(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);(p/'x').write_text('abc');write_tsv(p/'PACKAGE_MANIFEST.tsv',[{'path':'x','size_bytes':3,'sha256':sha_file(p/'x')}]);self.assertEqual(validate_manifest(p),1);(p/'x').write_text('abd')
   with self.assertRaises(ValueError):validate_manifest(p)
if __name__=='__main__':unittest.main()
