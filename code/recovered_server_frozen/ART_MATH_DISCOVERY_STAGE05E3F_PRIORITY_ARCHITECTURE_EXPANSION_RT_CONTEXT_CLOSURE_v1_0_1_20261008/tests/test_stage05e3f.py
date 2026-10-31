#!/usr/bin/env python3
from pathlib import Path
import tempfile, tarfile, csv, json, subprocess, hashlib, shutil

def wt(p,rows):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def make_fixture(td, with_rt_tar):
    td=Path(td); e=td/'e3e'; e.mkdir(); p=td/'project'; (p/'S/outputs').mkdir(parents=True)
    json.dump({'gate':'PASS_E3E_PRIORITY_FAMILY_LOCAL_HOMOLOGY_AND_REPEAT_COLOCALIZATION_CLOSURE','gate_pass':True},open(e/'STAGE05E3E_GATE.json','w'))
    hits=[]; col=[]; lock=[]
    famdefs=[('ARCH250_0010','F1','Q1','G1',True,1e-20,True,277,'UPSTREAM'),('ARCH250_0010','F1','Q1x','G2',False,1e-3,False,757,'UPSTREAM'),('ARCH250_0013','F2','Q2','G3',True,1e-30,True,85,'UPSTREAM'),('ARCH250_0013','F2','Q2x','G4',False,1e-20,True,110,'UPSTREAM'),('ARCH250_0033','F3','Q3','G5',True,1e-50,True,155,'DOWNSTREAM'),('ARCH250_0033','F3','Q3x','G6',False,1e-40,False,'','')]
    for fam,fc,rid,acc,wit,ev,close,gap,side in famdefs:
        hits.append({'family_id':fam,'fam50_cluster':fc,'target_name':rid,'sequence_id':acc,'full_evalue':ev,'full_score':100,'best_domain_evalue':ev,'best_domain_score':100,'target_aa_length':4,'global_alignment_identity':1 if wit else .5,'query_coverage_global':1,'target_coverage_global':1,'global_alignment_length':4,'exact_query_sequence':wit,'e3d_official_witness':wit,'target_aa_sha256':'x'})
        col.append({'family_id':fam,'target_name':rid,'sequence_id':acc,'full_evalue':ev,'e3d_official_witness':wit,'coordinate_resolved':True,'genome_resolved':True,'orf_start0':100,'orf_end0':112,'orf_strand':'+','repeat_core':'ACGT','repeat_core_length':4,'repeat_match_n_within_2kb_window':1 if gap!='' else 0,'nearest_repeat_found':gap!='','nearest_repeat_start0':50 if gap!='' else '','nearest_repeat_end0':54 if gap!='' else '','nearest_repeat_genomic_orientation':'+','nearest_repeat_relative_side':side,'nearest_repeat_gap_bp':gap,'exact_repeat_within_500bp':close,'exact_repeat_within_2000bp':gap!='','genome_fasta':'x','genome_sha256':'x'})
    for fam,fc in [('ARCH250_0010','F1'),('ARCH250_0013','F2'),('ARCH250_0033','F3')]: lock.append({'family_id':fam,'fam50_cluster':fc,'witness_n':1,'unique_witness_aa_n':1,'representative_aa_length':4,'representative_aa_sha256':'x','repeat_core_unique_n':1,'repeat_core_length':4,'repeat_core_sha256':'x','frozen_e3b_genomes':1,'frozen_e3b_taxonomy_tail_groups':1,'frozen_e3b_modal_side':'DOWNSTREAM' if fam.endswith('0033') else 'UPSTREAM','frozen_e3b_median_gap_bp':155 if fam.endswith('0033') else 100})
    wt(e/'STAGE05E3E_PHMMER_HITS.tsv',hits); wt(e/'STAGE05E3E_REPEAT_COLOCALIZATION.tsv',col); wt(e/'STAGE05E3E_PRIORITY_FAMILY_LOCK.tsv',lock)
    proteome=p/'S/outputs/FULL919_PREDICTED_PROTEINS.faa'; proteome.write_text(''.join(f'>{x[2]}\nAAAA\n' for x in famdefs))
    greg=p/'S/outputs/FULL919_GENE_REGISTRY.tsv'; wt(greg,[{'record_key':x[2],'sequence_id':x[3],'start0':100,'end0':112,'strand':'+'} for x in famdefs])
    json.dump({'full919_proteome_path':str(proteome),'full919_proteome_sha256':sha(proteome),'coordinate_registry':{'path':str(greg),'sha256':sha(greg)}},open(e/'STAGE05E3E_RESOURCE_STATUS.json','w'))
    if with_rt_tar:
        rtt=td/'rt.tsv'; wt(rtt,[{'sequence_id':'G4','start0':1000,'end0':1200,'record_key':'R1'}])
        tarpath=p/'ART_STAGE05C1_NEW919_BIOLOGY_RESULTS_20261006_092047.tar.gz'
        with tarfile.open(tarpath,'w:gz') as t: t.add(rtt,arcname='outputs/RT_LOCI_ALL919.tsv')
    etar=td/'e3e.tar.gz'
    with tarfile.open(etar,'w:gz') as t:
        for x in e.iterdir(): t.add(x,arcname=x.name)
    return p, etar

pkg=Path(__file__).resolve().parents[1]
# Exact-coordinate archive fallback still works.
with tempfile.TemporaryDirectory() as z:
    td=Path(z); p,etar=make_fixture(td,True)
    subprocess.run(['python',str(pkg/'scripts/stage05e3f.py'),'--e3e',str(etar),'--project-root',str(p),'--outdir',str(td/'out'),'--workdir',str(td/'work')],check=True,capture_output=True,text=True)
    g=json.load(open(td/'out/STAGE05E3F_GATE.json')); assert g['gate_pass'] and g['rt_context_resolved'] and g['rt_context_resolution_level']=='EXACT_COORDINATE_TABLE'
    s={r['family_id']:r for r in csv.DictReader(open(td/'out/STAGE05E3F_FAMILY_EXPANSION_SUMMARY.tsv'),delimiter='\t')}
    assert s['ARCH250_0013']['e3f_architecture_expansion_class']=='EXPANDED_SAME_ARCHITECTURE'
    assert s['ARCH250_0010']['e3f_architecture_expansion_class']=='WEAK_2KB_EXTENSION_ONLY'
    assert s['ARCH250_0033']['e3f_architecture_expansion_class']=='PROTEIN_FAMILY_BROAD_REPEAT_RESTRICTED'
# If exact RT coordinates are absent and homolog accessions have zero intersection
# with the frozen ten-RT-genome lock, accession-level closure is valid.
with tempfile.TemporaryDirectory() as z:
    td=Path(z); p,etar=make_fixture(td,False)
    subprocess.run(['python',str(pkg/'scripts/stage05e3f.py'),'--e3e',str(etar),'--project-root',str(p),'--outdir',str(td/'out'),'--workdir',str(td/'work')],check=True,capture_output=True,text=True)
    g=json.load(open(td/'out/STAGE05E3F_GATE.json')); assert g['gate_pass'] and g['rt_context_resolution_level']=='ACCESSION_LEVEL_ZERO_INTERSECTION' and g['frozen_rt_accession_lock_intersection_n']==0
    rr=list(csv.DictReader(open(td/'out/STAGE05E3F_RT_CONTEXT_SUMMARY.tsv'),delimiter='\t'))
    assert all(r['homologs_with_any_same_genome_stage05c1_rt_n']=='0' and r['within_10kb_of_stage05c1_rt_n']=='0' for r in rr)
print('[PASS] Stage05E3F exact-table + accession-lock-zero-intersection selftest')
