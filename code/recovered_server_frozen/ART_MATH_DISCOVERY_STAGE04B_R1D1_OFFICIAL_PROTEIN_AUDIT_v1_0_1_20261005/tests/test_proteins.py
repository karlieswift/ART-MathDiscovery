from pathlib import Path
import contextlib,copy,io,json,sys,tempfile,unittest,zipfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import audit_proteins as mod
from io_utils import digest,js,write_tsv,shatext

class ProteinTests(unittest.TestCase):
    def test_reverse_orientation_not_guessed_from_header(self):
        genome='AAAATGGCTTAAACC'
        nt=mod.rc(genome[3:12])
        parts=mod.parse_ranges({'accessionVersion':'X.1','range':[{'begin':'4','end':'12'}]},'X.1',len(genome))
        mapped,strand=mod.exact_cds_mapping(genome,nt,parts)
        self.assertEqual(strand,'-');self.assertEqual(mapped[0]['start0'],3)
    def test_join_parts_exact_nt(self):
        genome='ATGCCCTTAAAAGGG'
        parts=mod.parse_ranges({'accessionVersion':'X.1','range':[{'begin':'1','end':'3','orientation':'plus','order':0},{'begin':'7','end':'9','orientation':'plus','order':1}]},'X.1',len(genome))
        mapped,strand=mod.exact_cds_mapping(genome,'ATGTTA',parts)
        self.assertEqual(strand,'+');self.assertEqual(len(mapped),2)
    def test_endpoint_clipping_and_strict_overlap(self):
        p=[dict(start0=70,end0=80)]
        hit,d,lo,hi=mod.near_parts(p,90,100,10000,105)
        self.assertTrue(hit);self.assertEqual(hi,105);self.assertEqual(d,10)
    def test_wrong_nt_and_ambiguous_orientation_fail(self):
        with self.assertRaises(ValueError):mod.exact_cds_mapping('ATGC','ATGG',[dict(start0=0,end0=4,strand='?',order=None)])
        with self.assertRaises(ValueError):mod.exact_cds_mapping('ATAT','ATAT',[dict(start0=0,end0=4,strand='?',order=None)])
    def test_no_silent_stop_replacement(self):
        self.assertNotEqual(mod.protein_status('MA*AA'),'OFFICIAL_SEQUENCE')
        self.assertEqual(mod.protein_status('MAXAA'),'OFFICIAL_SEQUENCE')
    def test_exact_id_join_not_file_order(self):
        p={'P2.1:1-2':'MG','P1.1:1-2':'MA'}
        self.assertEqual(mod.select_record({'seqId':'P1.1:1-2'},p,'protein'),('P1.1:1-2','MA'))
        with self.assertRaises(ValueError):mod.select_record({'seqId':'P3.1'},p,'protein')
    def test_slippage_record_is_reported_not_silently_rebuilt(self):
        with self.assertRaises(ValueError):mod.parse_ranges({'accessionVersion':'X.1','range':[{'begin':1,'end':9,'ribosomalSlippage':-1}]},'X.1',10)
    def _run_full_synthetic(self,nested=False,orphan=False):
        # Full structural regression with artificial sequences, not real data.
        old_here,old_pol=mod.HERE,mod.POL
        try:
            with tempfile.TemporaryDirectory() as td:
                root=Path(td);mod.HERE=root/'tool';mod.HERE.mkdir();(mod.HERE/'inputs').mkdir();(mod.HERE/'outputs').mkdir()
                mod.POL=copy.deepcopy(old_pol)
                a04=root/'a04';r1a=root/'r1a';r1b=root/'r1b';r1c=root/'r1c'
                for p in (a04,r1a,r1b,r1c):(p/'outputs').mkdir(parents=True)
                (a04/'inputs').mkdir();(r1a/'work/zips').mkdir(parents=True)
                od=r1c/'outputs/stage04b_r1d_freeze';od.mkdir();of=r1c/'work/stage04b_r1d_freeze';of.mkdir(parents=True)
                candidates=[];stage=[];cache=[];old=[];old_fa=[]
                accs=[f'NC_{900000+i}.1' for i in range(173)]+['NC_138483.1']
                for i,acc in enumerate(accs):
                    num=0 if i==173 else 56 if i<45 else 55
                    genome='AC'+'ATGGCTTAA'*num+'A'*1000
                    gbytes=('>'+acc+'\n'+genome+'\n').encode()
                    frozen=a04/'inputs'/f'{acc}.frozen.fasta';frozen.write_bytes(gbytes)
                    stage.append(dict(sequence_id=acc,analysis_id=acc,frozen_fasta=str(frozen.relative_to(a04)),fasta_sha256=digest(frozen)))
                    for j in range(1 if i==173 else 3 if i<36 else 2):
                        candidates.append(dict(analysis_id=acc,sequence_id=acc,region_id=f'R{j+1:04d}',region_rank=str(j+1),length_bp=str(len(genome)),start0=str(600+j*10),end0=str(610+j*10)))
                    ann={'accession':acc};cdsfa=[];pfa=[];genes=[]
                    for j in range(num):
                        start=2+j*9;end=start+9;cid=f'{acc}:{start+1}-{end}';pid=f'P{i}_{j}.1:1-2'
                        cdsfa.append('>'+cid+'\nATGGCTTAA\n');pfa.append('>'+pid+'\nMA\n')
                        genes.append({'cds':[{'name':'hypothetical protein','nucleotide':{'accessionVersion':acc,'seqId':cid,'range':[{'begin':str(start+1),'end':str(end),'orientation':'plus'}]},'protein':{'seqId':pid,'accessionVersion':pid.split(':')[0]}}]})
                        pk=f'{acc}__{start}_{end}_+';aa='MX' if (i,j)==(0,0) else 'MA'
                        old.append(dict(protein_key=pk,sequence_id=acc,cds_start0=str(start),cds_end0=str(end),cds_strand='+',aa_length=len(aa),aa_sha256=shatext(aa)))
                        old_fa.append('>'+pk+'\n'+aa+'\n')
                    if nested and i==0:
                        children=[]
                        for child_index,(start,end,nt,aa) in enumerate([(2,5,'ATG','M'),(5,8,'GCT','A')]):
                            cid=f'{acc}:{start+1}-{end}';pid=f'CHILD{child_index}.1:1-1'
                            cdsfa.append('>'+cid+'\n'+nt+'\n');pfa.append('>'+pid+'\n'+aa+'\n')
                            children.append({'name':'synthetic processed peptide','nucleotide':{'seqId':cid,'accessionVersion':acc,'range':[{'begin':start+1,'end':end}]},'protein':{'seqId':pid,'accessionVersion':pid.split(':')[0]}})
                        genes[0]['cds'][0]['maturePeptide']=children[:1] if orphan else children
                    if genes:ann['genes']=genes
                    zpath=r1a/'work/zips'/f'{acc}.zip'
                    with zipfile.ZipFile(zpath,'w',compression=zipfile.ZIP_DEFLATED) as z:
                        z.writestr(mod.BASE+'genomic.fna',gbytes)
                        z.writestr(mod.BASE+'data_report.jsonl',json.dumps({'accession':acc,'virus':{'organismName':'SYNTHETIC'},'length':len(genome)})+'\n')
                        z.writestr(mod.BASE+'annotation_report.jsonl',json.dumps(ann)+'\n')
                        if num:z.writestr(mod.BASE+'cds.fna',''.join(cdsfa));z.writestr(mod.BASE+'protein.faa',''.join(pfa))
                    import hashlib
                    cache.append(dict(accession=acc,zip_sha256=digest(zpath),genomic_fna_sha256=hashlib.sha256(gbytes).hexdigest()))
                for i in range(66):stage.append(dict(sequence_id=f'EXTRA{i}',analysis_id=f'EXTRA{i}',frozen_fasta='unused',fasta_sha256='unused'))
                paths={
                    'stage04a_registry_sha256':mod.HERE/'inputs/STAGE04A_R1_DUAL_ARM_REGIONS.tsv',
                    'stage04a_universe_sha256':a04/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv',
                    'r1a_manifest_sha256':r1a/'outputs/STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv',
                    'r1d0_registry_sha256':od/'STAGE04B_R1D_NEIGHBORHOOD_PROTEIN_REGISTRY.tsv',
                    'r1d0_fasta_sha256':of/'STAGE04B_R1D_NEIGHBORHOOD_PROTEINS.faa'}
                for name,rows in [('stage04a_registry_sha256',candidates),('stage04a_universe_sha256',stage),('r1a_manifest_sha256',cache),('r1d0_registry_sha256',old)]:write_tsv(paths[name],rows,list(rows[0]))
                paths['r1d0_fasta_sha256'].write_text(''.join(old_fa))
                for k,p in paths.items():mod.POL[k]=digest(p)
                js(r1b/'outputs/STAGE04B_R1B_GATE.json',dict(gate_pass=True,n_exact_sequence_equivalent=174))
                with contextlib.redirect_stdout(io.StringIO()):result=mod.analyze(a04,r1a,r1b,r1c)
                self.assertEqual(result,2 if orphan else 0)
                gate=json.loads((mod.HERE/'outputs/GATE_STAGE04B_R1D1.json').read_text())
                self.assertEqual(gate['n_candidates_retained'],383);self.assertEqual(gate['n_old_translation_records'],9560)
                self.assertEqual(gate['n_official_neighborhood_proteins'],9560+(1 if orphan else 2 if nested else 0));self.assertEqual(gate['n_candidates_with_official_proteins'],382)
                self.assertEqual(gate['old_vs_official_counts']['AA_IDENTICAL'],9559)
                self.assertEqual(gate['old_vs_official_counts']['AA_DIFFERENT'],1)
                self.assertFalse(gate['negative_calls_authorized']);self.assertTrue(gate['orf_rescue_pending'])
                self.assertEqual(gate['official_protein_domain_scan_ready'],not orphan)
                if nested:
                    self.assertEqual(gate['n_neighborhood_mature_peptides'],1 if orphan else 2)
                    self.assertEqual(gate['n_mapping_issues'],2 if orphan else 0)
                    from io_utils import read_tsv
                    traversal=read_tsv(mod.HERE/'outputs/REPORT_NODE_TRAVERSAL.tsv')
                    children=[x for x in traversal if x['feature_type']=='MATURE_PEPTIDE']
                    self.assertTrue(all(x['parent_json_path']=='genes[0].cds[0]' for x in children))
                    if orphan:
                        issues=read_tsv(mod.HERE/'outputs/ANNOTATION_MAPPING_ISSUES.tsv')
                        self.assertEqual({x['kind'] for x in issues},{'UNLINKED_CDS_FASTA_RECORD','UNLINKED_PROTEIN_FASTA_RECORD'})
        finally:mod.HERE,mod.POL=old_here,old_pol

    def test_full_synthetic_174_genome_383_candidate_audit(self):
        self._run_full_synthetic()
    def test_full_synthetic_nested_peptides_recovered_without_dropping_candidates(self):
        self._run_full_synthetic(nested=True)
    def test_full_synthetic_orphan_remains_hold(self):
        self._run_full_synthetic(nested=True,orphan=True)
    def test_processed_peptide_tree_retains_parentage(self):
        a={'genes':[{'cds':[{'maturePeptide':[{'maturePeptide':[{}]}]}]}]}
        nodes=list(mod.iter_peptide_nodes(a))
        self.assertEqual([x[2] for x in nodes],['CDS','MATURE_PEPTIDE','MATURE_PEPTIDE'])
        self.assertEqual(nodes[2][4],'genes[0].cds[0].maturePeptide[0]')
    def test_unknown_strand_is_not_payload_plus(self):
        self.assertEqual(mod.declared_coding_strand({'nucleotide':{'range':[{'begin':1,'end':9}]}},{}),('?','NOT_DECLARED'))
    def test_gene_declared_minus_is_retained(self):
        gene={'nucleotide':{'range':[{'begin':1,'end':9,'orientation':'minus'}]}}
        self.assertEqual(mod.declared_coding_strand({'nucleotide':{'range':[{'begin':1,'end':9}]}},gene),('-','GENE_NUCLEOTIDE_METADATA'))
    def test_difference_patterns_are_not_cause_claims(self):
        self.assertEqual(mod.aa_difference_class('MA','MA'),'IDENTICAL')
        self.assertEqual(mod.aa_difference_class('VA','MA'),'FIRST_RESIDUE_ONLY')
        self.assertEqual(mod.aa_difference_class('MAA','MA'),'ONE_EXTRA_TERMINAL_RESIDUE_ONLY')
        self.assertEqual(mod.aa_difference_class('MXA','MR'),'NOT_FIRST_OR_TERMINAL_ONLY')
    def test_reference_inventory_sees_untraversed_fields_without_linking(self):
        refs=list(mod.fasta_reference_inventory({'other':[{'seqId':'X.1:1-9'}]}))
        self.assertEqual(refs,[{'json_path':'$.other[0].seqId','value':'X.1:1-9'}])
        self.assertEqual(list(mod.iter_peptide_nodes({'other':[{'seqId':'X.1:1-9'}]})),[])
    def test_malformed_processed_branch_fails_closed(self):
        with self.assertRaises(ValueError):list(mod.iter_peptide_nodes({'genes':[{'cds':[{'maturePeptide':{}}]}]}))

if __name__=='__main__':unittest.main()
