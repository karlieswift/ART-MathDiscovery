#!/usr/bin/env python3
import json,sys
from pathlib import Path
from common import *

def main():
    if len(sys.argv)!=6: raise SystemExit('usage finalize.py targets global art8 outdir lock')
    targ,gp,ap,outdir,lockp=sys.argv[1:]; out=Path(outdir)
    T=read_tsv(targ); G={r['family_rep_id']:r for r in read_tsv(gp)}; A={r['family_rep_id']:r for r in read_tsv(ap)}
    rows=[]
    for t in T:
        f=t['family_rep_id']; g=G[f]; a=A[f]
        public_art8=b(a['public_art8_rt_90pct_identity_80pct_query_coverage_compatible'])
        global_hi=b(g['rt_nr_high_identity_near_full_match']) or b(g['architecture_nt_near_full_match'])
        rows.append({
          'priority_rank':t['priority_rank'],'family_rep_id':f,'structure_tier':t['structure_tier'],'repeat_motif':t['repeat_motif'],'partner_phrog':t['partner_phrog'],
          'global_public_audit_state':g['global_public_audit_state'],'rt_nr_top_saccver':g['rt_nr_top_saccver'],'rt_nr_top_pident':g['rt_nr_top_pident'],'rt_nr_top_qcov':g['rt_nr_top_qcov'],'rt_nr_top_tcov':g['rt_nr_top_tcov'],'rt_nr_top_title':g['rt_nr_top_title'],
          'rt_nr_high_identity_near_full_match':g['rt_nr_high_identity_near_full_match'],'rt_nr_explicit_rt_annotation_top20':g['rt_nr_explicit_rt_annotation_top20'],'rt_nr_explicit_art_annotation_top50':g['rt_nr_explicit_art_annotation_top50'],
          'architecture_nt_near_full_match':g['architecture_nt_near_full_match'],'architecture_nt_top_saccver':g['architecture_nt_top_saccver'],'architecture_nt_top_pident':g['architecture_nt_top_pident'],'architecture_nt_top_qcov':g['architecture_nt_top_qcov'],
          'best_public_art8_accession':a['best_public_art8_accession'],'best_public_art8_name':a['best_public_art8_name'],'best_public_art8_art_id':a['best_public_art8_art_id'],'best_public_art8_rt_pident':a['best_rt_pident'],'best_public_art8_rt_qcov':a['best_rt_qcov'],
          'public_art8_rt_90pct_identity_80pct_query_coverage_compatible':a['public_art8_rt_90pct_identity_80pct_query_coverage_compatible'],'public_art8_architecture_near_full_match':a['public_art8_architecture_near_full_match'],
          'full_art95_crosswalk_status':'UNRESOLVED_FULL_ART95_MACHINE_READABLE_CATALOGUE_NOT_AVAILABLE',
          'absolute_global_novelty_claim_authorized':False,'anthropic_unreported_art_claim_authorized':False,
          'computational_novelty_closeout_status':('PUBLIC_ART8_90PCT_RT_SEQUENCE_COMPATIBLE' if public_art8 else ('PUBLIC_HIGH_ID_SEQUENCE_OR_LOCUS_PRESENT_OUTSIDE_ART8' if global_hi else 'NO_HIGH_ID_PUBLIC_MATCH_IN_NCBI_AUDIT_BUT_ART95_UNRESOLVED'))
        })
    fp=out/'F3B4E_FINAL_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_5.tsv'; write_tsv(fp,rows)
    gate={'stage':'Stage05F3B4E','gate':'PASS_F3B4E_GLOBAL_PUBLIC_NOVELTY_AUDIT_AND_PUBLIC_ART8_CROSSWALK_COMPLETE_WITH_FULL_ART95_BOUNDARY','gate_pass':True,
          'target_family_n':5,'public_ncbi_global_audit_complete':True,'public_named_art8_crosswalk_complete':True,'full_art95_crosswalk_complete':False,
          'absolute_global_novelty_claim_authorized':False,'anthropic_unreported_art_claim_authorized':False,
          'final_board_sha256':sha256_file(fp),'claim_boundary':'The two computational gaps are reduced to dated NCBI nr/nt evidence plus exact public ART8 crosswalk. Absolute novelty and absence from the unpublished/non-machine-readable ART95 remainder remain unresolvable without the full catalogue.'}
    (out/'STAGE05F3B4E_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2)); print('\n=== FINAL GLOBAL NOVELTY / ART CROSSWALK 5 ==='); print(fp.read_text().rstrip())

if __name__=='__main__': main()
