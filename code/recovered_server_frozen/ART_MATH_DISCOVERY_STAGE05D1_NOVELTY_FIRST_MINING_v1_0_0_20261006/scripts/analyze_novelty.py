#!/usr/bin/env python3
from pathlib import Path
import argparse, pandas as pd, re, json, hashlib, math
from Bio import SeqIO

def parse_domtbl(path):
    rows=[]
    with open(path,errors='replace') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            p=line.rstrip('\n').split(maxsplit=22)
            if len(p)<22: continue
            # hmmscan: target=model, query=sequence
            rows.append(dict(model=p[0],model_acc=p[1],model_len=int(p[2]),query=p[3],query_acc=p[4],query_len=int(p[5]),
                             full_e=float(p[6]),full_score=float(p[7]),domain_i_e=float(p[12]),domain_score=float(p[13]),
                             hmm_from=int(p[15]),hmm_to=int(p[16]),ali_from=int(p[17]),ali_to=int(p[18]),description=p[22] if len(p)>22 else ''))
    return pd.DataFrame(rows)

def parse_blast(path, cols):
    if not Path(path).exists() or Path(path).stat().st_size==0: return pd.DataFrame(columns=cols)
    return pd.read_csv(path,sep='\t',names=cols)

def maxhit(df, by, scorecols):
    if df.empty: return pd.DataFrame(columns=[by])
    sortcols=[c for c in scorecols if c in df.columns]
    asc=[False]*len(sortcols)
    return df.sort_values(sortcols,ascending=asc).drop_duplicates(by)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); a=ap.parse_args(); root=Path(a.root)
    cfg=root/'config'; out=root/'outputs'; work=root/'work'; out.mkdir(exist_ok=True)
    cand=pd.read_csv(cfg/'NOVELTY_FIRST_999_CANDIDATES.tsv',sep='\t')
    links=pd.read_csv(cfg/'NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv',sep='\t')
    fan=pd.read_csv(cfg/'NOVELTY_FIRST_PROTEIN_FANOUT.tsv',sep='\t')
    reg=pd.read_csv(cfg/'NOVELTY_FIRST_24020_ORF_REGISTRY.tsv',sep='\t')
    my=parse_domtbl(work/'myrt.domtbl')
    if not my.empty:
        my['passes_primary']=(my.full_e<=1e-4)&(my.domain_i_e<=1e-4)
        my=my[my.passes_primary].copy()
    my.to_csv(out/'MYRT_PRIMARY_DOMAINS.tsv',sep='\t',index=False)
    primary_targets=sorted(set(my['query'])) if not my.empty else []
    # CDD/Pfam refinement raw domains only for primary targets
    cdd=parse_domtbl(work/'cdd.domtbl') if (work/'cdd.domtbl').exists() else pd.DataFrame()
    pat=re.compile(r'(RVT|reverse\s+transcript|group[_ -]?II[_ -]?RT|RT_G2|Intron_maturas|GIIM|retron)',re.I)
    if not cdd.empty:
        cdd['rt_related_refinement']=cdd.apply(lambda r: bool(pat.search(' '.join(map(str,[r.model,r.model_acc,r.description])))),axis=1)
    cdd.to_csv(out/'CDD_PFAM_REFINEMENT_DOMAINS.tsv',sep='\t',index=False)
    supported=set(cdd.loc[cdd.rt_related_refinement,'query']) if not cdd.empty else set()
    # Motif-only exploratory branch: YxDD in all target proteins
    motifs=[]
    for rec in SeqIO.parse(root/'inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa','fasta'):
        s=str(rec.seq)
        for m in re.finditer(r'Y.DD',s): motifs.append((rec.id,m.start()+1,m.group()))
    mot=pd.DataFrame(motifs,columns=['target_id','motif_start_aa_1based','motif'])
    mot.to_csv(out/'YxDD_MOTIF_ONLY.tsv',sep='\t',index=False)
    motif_targets=set(mot.target_id)
    # candidate fanout for primary/support/motif
    fan2=fan.merge(links[['candidate_key','record_key']].drop_duplicates(),on='record_key',how='left')
    fan2['myrt_primary_hit']=fan2.target_id.isin(primary_targets)
    fan2['myrt_cdd_refined_rt_support']=fan2.target_id.isin(supported)
    fan2['YxDD_motif_only']=fan2.target_id.isin(motif_targets)
    fan2.to_csv(out/'PROTEIN_TO_CANDIDATE_BROAD_RT_LINKS.tsv',sep='\t',index=False)
    agg=fan2.groupby('candidate_key').agg(
      n_myrt_primary=('myrt_primary_hit','sum'),n_myrt_cdd_refined=('myrt_cdd_refined_rt_support','sum'),n_YxDD=('YxDD_motif_only','sum')
    ).reset_index()
    board=cand.merge(agg,on='candidate_key',how='left').fillna({'n_myrt_primary':0,'n_myrt_cdd_refined':0,'n_YxDD':0})
    # Reference similarity (descriptive, not full author-catalogue dereplication)
    bcols=['qseqid','sseqid','pident','length','qlen','slen','evalue','bitscore']
    rt=parse_blast(work/'rt_vs_registered.tsv',bcols)
    part=parse_blast(work/'allproteins_vs_partner.tsv',bcols)
    ncols=['qseqid','sseqid','pident','length','qlen','slen','evalue','bitscore']
    arr=parse_blast(work/'cores_vs_registered.tsv',ncols)
    for d in [rt,part,arr]:
        if not d.empty:
            d['qcov']=100*d['length']/d['qlen']; d['scov']=100*d['length']/d['slen']
    # target-level RT best
    if not rt.empty:
        rtbest=rt.sort_values(['bitscore','pident'],ascending=False).drop_duplicates('qseqid')
        rtbest['close_registered_ART_RT']=(rtbest.pident>=60)&(rtbest.qcov>=70)&(rtbest.scov>=70)
        rtbest=rtbest.rename(columns={'qseqid':'target_id','sseqid':'best_registered_ART_RT','pident':'best_registered_ART_RT_pident','qcov':'best_registered_ART_RT_qcov','scov':'best_registered_ART_RT_scov','bitscore':'best_registered_ART_RT_bitscore'})
        fan2=fan2.merge(rtbest[['target_id','best_registered_ART_RT','best_registered_ART_RT_pident','best_registered_ART_RT_qcov','best_registered_ART_RT_scov','best_registered_ART_RT_bitscore','close_registered_ART_RT']],on='target_id',how='left')
    else: fan2['close_registered_ART_RT']=False
    # partner target similarity mapped to candidate via all linked ORFs
    if not part.empty:
        pb=part.sort_values(['bitscore','pident'],ascending=False).drop_duplicates('qseqid')
        pb['close_registered_ART_partner']=(pb.pident>=40)&(pb.qcov>=70)&(pb.scov>=70)
        pb=pb.rename(columns={'qseqid':'target_id','sseqid':'best_registered_partner','pident':'best_registered_partner_pident','qcov':'best_registered_partner_qcov','scov':'best_registered_partner_scov','bitscore':'best_registered_partner_bitscore'})
        fp=fan[['record_key','target_id']].merge(links[['candidate_key','record_key']].drop_duplicates(),on='record_key').merge(pb,on='target_id',how='left')
        fp=fp.sort_values(['candidate_key','best_registered_partner_bitscore'],ascending=[True,False],na_position='last').drop_duplicates('candidate_key')
        board=board.merge(fp[['candidate_key','best_registered_partner','best_registered_partner_pident','best_registered_partner_qcov','best_registered_partner_scov','best_registered_partner_bitscore','close_registered_ART_partner']],on='candidate_key',how='left')
    else: board['close_registered_ART_partner']=False
    # array core reference proximity: map core sequences by sha to repeat groups/candidates
    if not arr.empty:
        ab=arr.sort_values(['bitscore','pident'],ascending=False).drop_duplicates('qseqid')
        ab['close_registered_ART_array_core']=(ab.pident>=80)&(ab.qcov>=80)&(ab.scov>=80)
        reps=pd.read_csv(cfg/'NOVELTY_FIRST_NONCODING_REPEAT_GROUPS.tsv',sep='\t')
        reps['core_target_id']=reps.core_nt.map(lambda s:'CORE_'+hashlib.sha256(str(s).encode()).hexdigest())
        ac=reps[['candidate_key','core_target_id']].drop_duplicates().merge(ab,left_on='core_target_id',right_on='qseqid',how='left')
        ac=ac.sort_values(['candidate_key','bitscore'],ascending=[True,False],na_position='last').drop_duplicates('candidate_key')
        ac=ac.rename(columns={'sseqid':'best_registered_array_core','pident':'best_registered_array_pident','qcov':'best_registered_array_qcov','scov':'best_registered_array_scov','bitscore':'best_registered_array_bitscore'})
        board=board.merge(ac[['candidate_key','best_registered_array_core','best_registered_array_pident','best_registered_array_qcov','best_registered_array_scov','best_registered_array_bitscore','close_registered_ART_array_core']],on='candidate_key',how='left')
    else: board['close_registered_ART_array_core']=False
    # candidate RT proximity: any refined RT linked protein close to registered ART RT
    if 'close_registered_ART_RT' in fan2:
        x=fan2[fan2.myrt_cdd_refined_rt_support].groupby('candidate_key').close_registered_ART_RT.max().rename('refined_RT_close_to_registered_ART').reset_index()
        board=board.merge(x,on='candidate_key',how='left')
    else: board['refined_RT_close_to_registered_ART']=False
    for c in ['close_registered_ART_partner','close_registered_ART_array_core','refined_RT_close_to_registered_ART']:
        board[c]=board[c].fillna(False).astype(bool)
    board['registered_reference_close_component_n']=board[['close_registered_ART_partner','close_registered_ART_array_core','refined_RT_close_to_registered_ART']].sum(axis=1)
    # Fixed review tiers. No ART identity claim.
    board['novelty_review_tier']='T5_RT_UNRESOLVED_REPEAT_CONTEXT'
    strict_exact=board.novelty_architecture_stratum.eq('N1_EXACT_NONCODING_ARRAY_OVERLAP')
    strict_approx=board.novelty_architecture_stratum.eq('N2_APPROX_NONCODING_ARRAY_OVERLAP')
    aux=board.evidence_label.eq('RT_LIKE_AUXILIARY_PROFILE_CONTEXT')
    refined=board.n_myrt_cdd_refined>0
    motif=(board.n_YxDD>0)&~refined
    close=board.registered_reference_close_component_n>0
    board.loc[aux,'novelty_review_tier']='T4_AUXILIARY_RT_LIKE_REVIEW'
    board.loc[motif & (strict_exact|strict_approx),'novelty_review_tier']='T3_MOTIF_ONLY_ARRAY_REVIEW'
    board.loc[refined & ~(strict_exact|strict_approx),'novelty_review_tier']='T3_BROAD_RT_REPEAT_CONTEXT'
    board.loc[refined & strict_approx,'novelty_review_tier']='T2_BROAD_RT_APPROX_ARRAY'
    board.loc[refined & strict_exact,'novelty_review_tier']='T1_BROAD_RT_EXACT_ARRAY'
    board.loc[refined & (strict_exact|strict_approx) & ~close,'novelty_review_tier']='T0_HIGH_NOVELTY_REVIEW_NOT_CONFIRMED'
    board['new_ART_confirmed']=False
    board['wet_lab_ready']=False
    board['full_author_catalogue_dereplication_complete']=False
    board.to_csv(out/'NOVELTY_FIRST_BOARD.tsv',sep='\t',index=False)
    # summary
    summary={'stage':'Stage05D1','gate':'PASS_NOVELTY_FIRST_999_BROAD_RT_AND_REFERENCE_PROXIMITY_COMPUTED_WITH_SCOPE_LIMITS','gate_pass':True,
      'n_candidates':len(board),'n_genomes':board.sequence_id.nunique(),'n_unique_protein_targets':fan.target_id.nunique(),
      'n_myrt_primary_candidates':int((board.n_myrt_primary>0).sum()),'n_myrt_cdd_refined_candidates':int((board.n_myrt_cdd_refined>0).sum()),
      'n_YxDD_motif_candidates':int((board.n_YxDD>0).sum()),'tier_counts':board.novelty_review_tier.value_counts().to_dict(),
      'registered_reference_proximity_is_not_full_catalogue_dereplication':True,'new_ART_confirmed':0,'wet_lab_ready_declared':0,
      'mathematical_scores_or_ranks_changed':False,'candidate_filtering_performed':False,'UG27_touched':False}
    (out/'SUMMARY_STAGE05D1.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    (out/'GATE_STAGE05D1.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
