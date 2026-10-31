#!/usr/bin/env python3
import argparse, json, hashlib, os, re
from pathlib import Path
import pandas as pd
import numpy as np

EXPECTED_SHA='31c70ee1c7769baec58da11dfcc40fc3ec6ab7c17d6608391f737bc404041bca'
FLAGCOLS=['strict_rt_keyword_10kb_n','retron_keyword_10kb_n','dgr_keyword_10kb_n','crispr_keyword_10kb_n','groupII_keyword_10kb_n','telomerase_keyword_10kb_n','retroelement_keyword_10kb_n']

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def lineage_class(r):
    if int(r['genome_n'])>=5 and int(r['unique_taxonomy_tail3_n'])>=3 and float(r['modal_taxonomy_tail3_fraction'])<=0.75:
        return 'BROAD_MULTI_LINEAGE'
    if int(r['genome_n'])>=3 and int(r['unique_taxonomy_tail3_n'])>=2:
        return 'MULTI_LINEAGE'
    return 'NARROW_LINEAGE_REPLICATED'

def review_stratum(r):
    broad=r['lineage_breadth_class']
    sup=float(r['official_same_strand_high_overlap_or_near_exact_witness_fraction'])
    side=float(r['modal_side_fraction_E2'])
    kw=int(r['known_system_keyword_candidate_n'])
    if broad=='BROAD_MULTI_LINEAGE' and sup>=0.75 and side>=0.75 and kw==0:
        return 'E3B1_BROAD_COHERENT_REPLICATED'
    if broad=='MULTI_LINEAGE' and sup>=0.75 and side>=0.75 and kw==0:
        return 'E3B2_MULTI_LINEAGE_COHERENT_REPLICATED'
    if broad in ('BROAD_MULTI_LINEAGE','MULTI_LINEAGE'):
        return 'E3B3_MULTI_LINEAGE_HETEROGENEOUS_RETAINED'
    if broad=='NARROW_LINEAGE_REPLICATED' and sup>=0.50 and side>=0.75 and kw==0:
        return 'E3B4_NARROW_COHERENT_REPLICATED'
    return 'E3B5_SPECIAL_CONFLICT_REVIEW'

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--e3a-root', required=True)
    ap.add_argument('--outdir', required=True)
    a=ap.parse_args()
    root=Path(a.e3a_root); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    fam=pd.read_csv(root/'outputs/STAGE05E3A_FAMILY_BIOLOGICAL_CONTEXT.tsv',sep='\t')
    e2fam=pd.read_csv(root/'frozen/STAGE05E3A_E2A_FAMILIES.tsv',sep='\t')
    cand=pd.read_csv(root/'outputs/STAGE05E3A_CANDIDATE_OFFICIAL_CONTEXT.tsv',sep='\t')
    neigh=pd.read_csv(root/'outputs/STAGE05E3A_NEIGHBORHOOD_OFFICIAL_CDS.tsv',sep='\t')
    meta=pd.read_csv(root/'outputs/STAGE05E3A_SEQUENCE_METADATA.tsv',sep='\t')
    wit=pd.read_csv(root/'frozen/STAGE05E3A_E2A_WITNESS_PAIRS.tsv',sep='\t')
    gate=json.load(open(root/'outputs/GATE_STAGE05E3A.json'))
    if not gate.get('gate_pass'): raise SystemExit('E3A gate not PASS')
    if (len(fam),len(cand),cand.sequence_id.nunique(),len(wit))!=(12,119,61,80):
        raise SystemExit(f'Frozen count mismatch: {(len(fam),len(cand),cand.sequence_id.nunique(),len(wit))}')

    # family-level host diversity proxy
    fc=cand[['family_id','sequence_id']].drop_duplicates().merge(meta,on='sequence_id',how='left')
    fc['host_genus']=fc['host'].fillna('').astype(str).str.split().str[0].replace('',np.nan)
    hostrows=[]
    for fid,g in fc.groupby('family_id'):
        hg=g['host_genus'].fillna('NA')
        hv=hg.value_counts()
        hostrows.append({'family_id':fid,'unique_host_genus_n':int(hg.nunique()),'modal_host_genus':hv.index[0],'modal_host_genus_fraction':float(hv.iloc[0]/len(g))})
    hostdf=pd.DataFrame(hostrows)

    # candidate keyword burden per family and exact categories
    kwrows=[]
    for fid,g in cand.groupby('family_id'):
        flags=(g[FLAGCOLS].fillna(0).astype(int)>0)
        kwrows.append({'family_id':fid,'known_system_keyword_candidate_n':int(flags.any(axis=1).sum())})
    kwdf=pd.DataFrame(kwrows)

    # matched official witness protein length + annotation
    w=cand[(cand.has_exact_primary_joint_witness.astype(int)==1)&cand.witness_official_match_class.fillna('').str.contains('SAME_STRAND')].copy()
    nmap=neigh[['candidate_key','protein_id','cds_start0','cds_end0','product']].drop_duplicates()
    wm=w.merge(nmap,left_on=['candidate_key','witness_official_protein_id'],right_on=['candidate_key','protein_id'],how='left',suffixes=('','_neigh'))
    wm['official_witness_aa_len']=(wm['cds_end0']-wm['cds_start0'])/3.0
    prod=wm['witness_official_product'].fillna('').str.lower()
    wm['official_witness_hypothetical_or_uncharacterized']=prod.str.contains('hypothetical|uncharacterized|unknown').astype(int)
    wm_out=wm[['family_id','candidate_key','sequence_id','organism','taxonomy_tail3','host','witness_official_protein_id','witness_official_locus_tag','witness_official_product','official_witness_aa_len','official_witness_hypothetical_or_uncharacterized','witness_official_overlap_fraction']].copy()
    wm_out.to_csv(out/'STAGE05E3B_OFFICIAL_WITNESS_PROTEINS.tsv',sep='\t',index=False)

    mrows=[]
    for fid,g in wm.groupby('family_id'):
        mrows.append({'family_id':fid,'official_supported_witness_n':len(g),'official_witness_median_aa':float(g.official_witness_aa_len.median()),'official_witness_min_aa':float(g.official_witness_aa_len.min()),'official_witness_max_aa':float(g.official_witness_aa_len.max()),'official_witness_hypothetical_fraction':float(g.official_witness_hypothetical_or_uncharacterized.mean())})
    mdf=pd.DataFrame(mrows)

    board=fam.merge(e2fam[['family_id','joint_exact_primary_signature_edge_fraction','joint_any_signature_edge_fraction','family_exclusive_exact_primary_core_n','family_exclusive_compact_FAM50_n']],on='family_id',how='left').merge(hostdf,on='family_id',how='left').merge(kwdf,on='family_id',how='left').merge(mdf,on='family_id',how='left')
    board['lineage_breadth_class']=board.apply(lineage_class,axis=1)
    board['stage05e3b_review_stratum']=board.apply(review_stratum,axis=1)
    order={'E3B1_BROAD_COHERENT_REPLICATED':1,'E3B2_MULTI_LINEAGE_COHERENT_REPLICATED':2,'E3B3_MULTI_LINEAGE_HETEROGENEOUS_RETAINED':3,'E3B4_NARROW_COHERENT_REPLICATED':4,'E3B5_SPECIAL_CONFLICT_REVIEW':5}
    board['stage05e3b_review_order_group']=board.stage05e3b_review_stratum.map(order)
    board=board.sort_values(['stage05e3b_review_order_group','review_order','family_id'])
    board.to_csv(out/'STAGE05E3B_FAMILY_BOARD.tsv',sep='\t',index=False)
    board.to_csv(out/'STAGE05E3B_REVIEW_QUEUE.tsv',sep='\t',index=False)

    # candidate-specific known system context (specific labels; CRISPR is defense context, not retro)
    k=neigh[(neigh[['strict_rt_keyword','flag_retron','flag_dgr','flag_crispr','flag_groupII','flag_telomerase','flag_retroelement']].fillna(0).astype(int)>0).any(axis=1)].copy()
    if len(k):
        def labels(r):
            xs=[]
            if int(r.strict_rt_keyword): xs.append('STRICT_RT_ANNOTATION')
            if int(r.flag_retron): xs.append('RETRON_CONTEXT')
            if int(r.flag_dgr): xs.append('DGR_CONTEXT')
            if int(r.flag_crispr): xs.append('CRISPR_DEFENSE_CONTEXT')
            if int(r.flag_groupII): xs.append('GROUPII_INTRON_CONTEXT')
            if int(r.flag_telomerase): xs.append('TELOMERASE_CONTEXT')
            if int(r.flag_retroelement): xs.append('RETROELEMENT_CONTEXT')
            return ';'.join(xs)
        k['context_labels']=k.apply(labels,axis=1)
    k.to_csv(out/'STAGE05E3B_KNOWN_SYSTEM_CONTEXT.tsv',sep='\t',index=False)

    # candidate flags all retained
    cf=cand[['family_id','candidate_key','sequence_id','region_rank','novelty_architecture_stratum','has_exact_primary_joint_witness','witness_official_match_class','witness_official_product','organism','taxonomy_tail3','host']].copy()
    cf=cf.merge(board[['family_id','lineage_breadth_class','stage05e3b_review_stratum']],on='family_id',how='left')
    cf.to_csv(out/'STAGE05E3B_CANDIDATE_FLAGS.tsv',sep='\t',index=False)

    # summary
    strata=board.stage05e3b_review_stratum.value_counts().to_dict()
    lineages=board.lineage_breadth_class.value_counts().to_dict()
    total_w=len(wm)
    hyp_n=int(wm.official_witness_hypothetical_or_uncharacterized.sum())
    le100=int((wm.official_witness_aa_len<=100).sum())
    summary={
      'stage':'Stage05E3B','n_families':12,'n_candidates':119,'n_sequence_ids':61,'n_exact_primary_joint_witness_candidates':80,
      'n_official_supported_witness_candidates':total_w,
      'official_supported_witness_fraction_of_exact_primary_witness':total_w/80.0,
      'n_official_supported_witness_hypothetical_or_uncharacterized':hyp_n,
      'official_supported_witness_hypothetical_fraction':hyp_n/total_w,
      'official_supported_witness_median_aa':float(wm.official_witness_aa_len.median()),
      'n_official_supported_witness_le_100aa':le100,
      'official_supported_witness_le_100aa_fraction':le100/total_w,
      'lineage_breadth_counts':lineages,
      'review_stratum_counts':strata,
      'n_families_with_strict_rt_keyword_10kb':int((board.strict_rt_keyword_candidate_n_10kb>0).sum()),
      'n_families_with_any_known_system_keyword_10kb':int((board.known_system_keyword_candidate_n>0).sum()),
      'candidate_filtering_performed':False,'family_membership_changed':False,'selection_uses_structure_results':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,
      'gate':'PASS_E2A_FAMILY_LINEAGE_DEREPLICATION_AND_CONTEXT_SYNTHESIS_FROZEN',
      'gate_pass':True,
      'claim_boundary':'E3B organizes all 12 frozen E2A families by lineage breadth, topology coherence, official CDS support, and explicit known-system context. It does not establish ART identity or novelty.'
    }
    # Regression expectations protect implementation semantics
    expected_strata={'E3B1_BROAD_COHERENT_REPLICATED':1,'E3B2_MULTI_LINEAGE_COHERENT_REPLICATED':2,'E3B3_MULTI_LINEAGE_HETEROGENEOUS_RETAINED':5,'E3B4_NARROW_COHERENT_REPLICATED':2,'E3B5_SPECIAL_CONFLICT_REVIEW':2}
    expected_lineages={'BROAD_MULTI_LINEAGE':4,'MULTI_LINEAGE':4,'NARROW_LINEAGE_REPLICATED':4}
    checks=[strata==expected_strata,lineages==expected_lineages,total_w==64,hyp_n==63,le100==52,int((board.strict_rt_keyword_candidate_n_10kb>0).sum())==0,int((board.known_system_keyword_candidate_n>0).sum())==1]
    if not all(checks):
        summary['gate']='FAIL_STAGE05E3B_REGRESSION_MISMATCH'; summary['gate_pass']=False
    json.dump(summary,open(out/'STAGE05E3B_SUMMARY.json','w'),indent=2)
    json.dump(summary,open(out/'GATE_STAGE05E3B.json','w'),indent=2)
    # compact discoveries table, descriptive not selection
    disc=board[['stage05e3b_review_stratum','review_order','family_id','candidate_n','genome_n','unique_taxonomy_tail3_n','unique_host_genus_n','signature_witness_candidate_n','official_same_strand_high_overlap_or_near_exact_witness_n','official_same_strand_high_overlap_or_near_exact_witness_fraction','official_witness_median_aa','official_witness_hypothetical_fraction','median_repeat_orf_gap_bp_E2','modal_transcriptional_repeat_side_E2','modal_side_fraction_E2','strict_rt_keyword_candidate_n_10kb','known_system_keyword_candidate_n']].copy()
    disc.to_csv(out/'STAGE05E3B_DISCOVERY_SUMMARY.tsv',sep='\t',index=False)
    # hashes
    hs=[]
    for p in sorted(out.iterdir()):
        if p.is_file() and p.name!='STAGE05E3B_SHA256.tsv': hs.append((p.name,sha256_file(p)))
    pd.DataFrame(hs,columns=['file','sha256']).to_csv(out/'STAGE05E3B_SHA256.tsv',sep='\t',index=False)
    print(json.dumps(summary,indent=2))
    if not summary['gate_pass']: raise SystemExit(2)

if __name__=='__main__': main()
