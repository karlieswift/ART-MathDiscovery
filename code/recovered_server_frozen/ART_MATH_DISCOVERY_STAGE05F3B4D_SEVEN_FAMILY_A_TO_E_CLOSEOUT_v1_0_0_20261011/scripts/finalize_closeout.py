#!/usr/bin/env python3
import csv,json,sys,hashlib
from pathlib import Path

def boolv(x):return str(x).lower()=='true'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):h.update(b)
    return h.hexdigest()
def rank_tier(t):
    return {'STRONG':0,'SUGGESTIVE':1,'WEAK':2,'NO_HIT':3,'TECHNICAL_FAILURE':4,'FOLDSEEK_TECHNICAL_FAILURE':4}.get(t,5)
def readmap(p,key):
    return {r[key]:r for r in csv.DictReader(open(p,encoding='utf-8'),delimiter='\t')}

def main():
    if len(sys.argv)!=13:
        raise SystemExit('usage: finalize_closeout.py TARGET7 A_SUM B_SUM C_SUM D_SUM E_SUM A_GATE B_GATE C_GATE D_GATE E_GATE OUTDIR')
    target_p,a_p,b_p,c_p,d_p,e_p,ag,bg,cg,dg,eg,outdir=sys.argv[1:]
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    gates=[]
    for p in [ag,bg,cg,dg,eg]:
        g=json.load(open(p));
        if not g.get('gate_pass'):raise SystemExit(f'FAIL upstream gate {p}')
        gates.append(g)
    targets=list(csv.DictReader(open(target_p,encoding='utf-8'),delimiter='\t'))
    if len(targets)!=7:raise SystemExit('FAIL target 7')
    A=readmap(a_p,'family_rep_id');B=readmap(b_p,'family_rep_id');C=readmap(c_p,'family_rep_id');E=readmap(e_p,'family_rep_id')
    D=readmap(d_p,'family_rep_id')
    rows=[]
    for t in targets:
        fam=t['family_rep_id']
        if boolv(t['preexisting_structure_tested']):
            st=t['preexisting_structure_tier']; scid=t['preexisting_structure_candidate_id']; source='PREEXISTING_F3B3E_R1'
        else:
            dr=D.get(fam,{})
            st=dr.get('r5c_tier','TECHNICAL_FAILURE');scid=dr.get('candidate_id','');source='F3B4D_D_STRUCTURE_EXTENSION'
        if st in {'STRONG','SUGGESTIVE'}:
            comp='RT_FOLD_SUPPORTED_STRICT_ARRAY_RECURRENT_PARTNER_CANDIDATE'
        elif st in {'WEAK','NO_HIT'}:
            comp='STRICT_ARRAY_RECURRENT_PARTNER_WITHOUT_FROZEN_RT_FOLD_SUPPORT'
        else:
            comp='STRUCTURE_UNRESOLVED_STRICT_ARRAY_RECURRENT_PARTNER_CANDIDATE'
        novelty=E[fam]['novelty_status']
        z={
          'family_rep_id':fam,
          'structure_tier':st,'structure_candidate_id':scid,'structure_evidence_source':source,
          'strict_array_locus_n':t['strict_array_locus_n'],'strict_recurrent_partner_accession_n':t['strict_array_dominant_partner_cluster_accession_n'],'strict_recurrent_partner_median_gap_bp':t['strict_array_dominant_partner_cluster_median_gap_bp'],
          'architecture_core_locus_n':A[fam]['core_architecture_locus_n'],'architecture_near_context_cluster_n':A[fam]['near_context_cluster_n'],'architecture_exact_context_cluster_n':A[fam]['exact_architecture_context_cluster_n'],
          'dominant_repeat_motif_neighborhood':B[fam]['dominant_motif_neighborhood_members'],'modal_exact_motif_core':B[fam]['modal_exact_motif_core'],'spacer_bp_median':B[fam]['spacer_bp_median'],'copy_n_median':B[fam]['copy_n_median'],
          'partner_exact_aa_sequence_n':C[fam]['exact_partner_aa_sequence_n'],'partner_pairwise_global_identity_min':C[fam].get('partner_pairwise_global_identity_min',''),'partner_pairwise_global_identity_median':C[fam].get('partner_pairwise_global_identity_median',''),'partner_aa_length_median':C[fam]['partner_aa_length_median'],'partner_same_strand_n':C[fam]['same_strand_n'],'partner_opposite_strand_n':C[fam]['opposite_strand_n'],
          'partner_pfam_top_target':C[fam].get('pfam_top_target',''),'partner_phrog_top_target':C[fam].get('phrog_top_target',''),'partner_vog_top_target':C[fam].get('vog_top_target',''),
          'local_known_exact_full_length_match':E[fam]['local_exact_full_length_match'],'local_known_high_identity_near_full_match':E[fam]['local_high_identity_near_full_match'],'local_top_reference_id':E[fam]['local_top_reference_id'],'local_top_reference_explicit_rt_keyword':E[fam]['local_top_reference_explicit_rt_keyword'],'local_top_reference_public_art_named_subset_keyword':E[fam]['local_top_reference_public_art_named_subset_keyword'],
          'novelty_status':novelty,'anthropic_unreported_art_claim_authorized':False,
          'computational_art_like_status':comp,
          'biochemical_rt_activity_confirmed':False,'dedicated_partner_function_confirmed':False,
        }
        rows.append(z)
    rows.sort(key=lambda r:(rank_tier(r['structure_tier']),-int(r['strict_recurrent_partner_accession_n']),-int(r['architecture_near_context_cluster_n']),r['family_rep_id']))
    for i,r in enumerate(rows,1):r['computational_priority_rank']=i
    cols=['computational_priority_rank']+[c for c in rows[0] if c!='computational_priority_rank']
    board=out/'F3B4D_FINAL_SEVEN_FAMILY_COMPUTATIONAL_EVIDENCE_BOARD.tsv'
    with board.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    supported=sum(r['structure_tier'] in {'STRONG','SUGGESTIVE'} for r in rows)
    gate={
      'stage':'Stage05F3B4D','gate':'PASS_F3B4D_SEVEN_FAMILY_A_TO_E_COMPUTATIONAL_CLOSEOUT_COMPLETE','gate_pass':True,
      'target_family_n':7,'r5c_strong_or_suggestive_family_n':supported,'r5c_nonpositive_resolved_family_n':sum(r['structure_tier'] in {'WEAK','NO_HIT'} for r in rows),'structure_technical_unresolved_family_n':sum('TECHNICAL' in r['structure_tier'] for r in rows),
      'all_target_families_have_strict_array_and_recurrent_partner_by_definition':True,
      'A_near_context_derep_complete':True,'B_repeat_core_conservation_complete':True,'C_partner_conservation_complete':True,'D_structure_extension_complete':True,'E_known_novelty_crosswalk_complete_with_boundary':True,
      'biochemical_rt_activity_confirmed_family_n':0,'dedicated_partner_function_confirmed_family_n':0,
      'anthropic_unreported_art_claim_authorized':False,
      'evidence_board_sha256':sha(board),
      'claim_boundary':'Computational A-E closeout of seven strict-array + recurrent-nearest-partner RVT families. A family with STRONG/SUGGESTIVE R5C may be called a high-priority computational ART-like candidate system, but not a confirmed ART, not a biochemically confirmed RT system, and not Anthropic-unreported without complete global/ART95 novelty closure.'}
    (out/'STAGE05F3B4D_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2));print('\n=== FINAL SEVEN-FAMILY COMPUTATIONAL EVIDENCE BOARD ===');print('\t'.join(cols))
    for r in rows:print('\t'.join(str(r[c]) for c in cols))
    print('PASS_STAGE05F3B4D_SEVEN_FAMILY_A_TO_E_CLOSEOUT_COMPLETE')

if __name__=='__main__':main()
