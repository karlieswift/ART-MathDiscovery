"""Join physical and sequence evidence back to every frozen candidate."""
from __future__ import annotations
from collections import defaultdict,Counter
from common import *


def integrate(root,fan,domains,repeats,windows,all_orfs):
    root=Path(root);candidates=read_tsv(root/'inputs/candidates.tsv');bytarget=defaultdict(set)
    for d in domains:
        if d['evidence_tier']!='REPORTED_ONLY':bytarget[d['target_id']].add(d['evidence_tier'])
    bycand=defaultdict(list);assessed=[]
    for r in fan:
        tiers=bytarget[r['target_id']];p=dict(r)
        p.update({'rt_core_GA':int('RT_CORE_GA' in tiers),'rt_core_sensitivity_only':int('RT_CORE_GA' not in tiers and 'RT_CORE_SENSITIVITY' in tiers),'rt_like_auxiliary_GA':int('RT_LIKE_AUXILIARY_GA' in tiers),'rt_profile_status':';'.join(sorted(tiers)) if tiers else 'NO_HIT_IN_THIS_PROFILE_PANEL','RT_activity_proven':False})
        assessed.append(p)
        for ck in r['candidate_keys'].split(';'):
            if ck:bycand[ck].append(p)
    write_tsv(root/'outputs/ALL_PROTEIN_RT_ASSESSMENTS.tsv',assessed,list(assessed[0]))
    # Already-known textual RT case: diagnostic only, never a threshold-tuning target.
    write_tsv(root/'outputs/KNOWN_RT_TEXT_CASE_SCAN_STATUS.tsv',[r for r in assessed if r['official_protein_accession']=='YP_010110885.1'],list(assessed[0]))
    R=defaultdict(list)
    for r in repeats:R[r['candidate_key']].append(r)
    W={r['candidate_key']:r for r in windows};features=defaultdict(list)
    for r in read_tsv(root/'inputs/all_official_features.tsv'):
        if r['feature_type']=='CDS':features[r['sequence_id']].append(r)
    orfby=defaultdict(list)
    for r in all_orfs:orfby[r['sequence_id']].append(r)
    evidence=[];pairs=[];partners=[];seenpartners=set();tier_counts=Counter()
    for c in candidates:
        ck=c['sequence_id']+'::'+c['region_id'];proteins=bycand[ck];ww=W[ck];rr=R[ck]
        rt=[p for p in proteins if p['rt_core_GA']]
        weak=[p for p in proteins if p['rt_core_sensitivity_only']]
        high=[r for r in rr if r['high_complexity_evidence']]
        relnc=[r for r in high if r['annotation_relative_repeat_context']=='FULL_SPAN_OUTSIDE_ANNOTATED_CDS']
        stringent=[r for r in relnc if r['predicted_orf_assessed'] and r['predicted_cds_overlap_bp_span']==0]
        if rt and stringent:tag='RT_CORE_GA_AND_REPEAT_SPAN_OUTSIDE_OFFICIAL_AND_PREDICTED_CDS'
        elif rt and relnc:tag='RT_CORE_GA_AND_ANNOTATION_RELATIVE_REPEAT_BUT_PREDICTED_ORF_OVERLAP'
        elif rt:tag='RT_CORE_GA_WITHOUT_STRICT_NONCODING_REPEAT_EVIDENCE'
        elif weak:tag='RT_CORE_SENSITIVITY_ONLY'
        elif any(p['rt_like_auxiliary_GA'] for p in proteins):tag='RT_LIKE_AUXILIARY_ONLY'
        elif any('RT_LIKE_AUXILIARY_SENSITIVITY' in p['rt_profile_status'] for p in proteins):tag='RT_LIKE_AUXILIARY_SENSITIVITY_ONLY'
        else:tag='NO_RT_PROFILE_HIT_AT_FROZEN_TIERS'
        tier_counts[tag]+=1
        # The tier is a descriptive label, not a candidate-survival decision.
        row=dict(c);row.update({'candidate_key':ck,'n_official_proteins_scanned':sum(p['source_branch']=='OFFICIAL' for p in proteins),'n_predicted_proteins_scanned':sum(p['source_branch']!='OFFICIAL' for p in proteins),'n_official_RT_core_GA':sum(p['source_branch']=='OFFICIAL' for p in rt),'n_genome_only_rescue_RT_core_GA':sum(p['source_branch']=='PRODIGAL_RESCUE' for p in rt),'n_annotation_gap_sensitivity_ORF_RT_core_GA':sum(p['source_branch']=='PRODIGAL_SENSITIVITY' for p in rt),'n_RT_core_GA_unique_AA':len({p['target_id'] for p in rt}),'n_RT_core_sensitivity_only_records':len(weak),'n_exact_core_repeat_groups':ww['n_exact_core_repeat_groups'],'n_high_complexity_repeat_groups':ww['n_high_complexity_repeat_groups'],'n_high_complexity_full_span_outside_official_CDS':ww['n_high_complexity_full_span_outside_official_cds'],'n_high_complexity_full_span_outside_official_and_predicted_CDS':ww['n_high_complexity_full_span_outside_official_and_predicted_cds'],'recurrence_window_official_CDS_overlap_bp':ww['recurrence_official_cds_overlap_bp'],'recurrence_window_predicted_CDS_overlap_bp':ww['recurrence_predicted_cds_overlap_bp'],'evidence_label':tag,'official_annotation_gap':int(not ww['official_annotation_available']),'RNA_expression_status':'NOT_ASSESSED','ART_identity_status':'NOT_ASSESSED','novelty_status':'NOT_ASSESSED','negative_call_authorized':False})
        evidence.append(row)
        for p in rt:
            pp=json.loads(p['parts_json'])
            for r in rr:
                a=(r['start0'],r['end0'])
                pairs.append({'candidate_key':ck,'sequence_id':c['sequence_id'],'rt_record_key':p['record_key'],'rt_source_branch':p['source_branch'],'rt_official_accession':p['official_protein_accession'],'rt_feature_type':p['feature_type'],'rt_parent_json_path':p['parent_json_path'],'repeat_group_id':r['repeat_group_id'],'RT_to_repeat_span_distance_bp':min(distance(x,a) for x in pp),'RT_overlaps_repeat_span':int(any(overlaps(x,a) for x in pp)),'repeat_high_complexity':r['high_complexity_evidence'],'repeat_context':r['annotation_relative_repeat_context'],'repeat_predicted_CDS_overlap_bp':r['predicted_cds_overlap_bp_span'],'repeat_overlaps_original_candidate':r['overlaps_frozen_candidate'],'interpretation':'same fixed neighborhood; not a functional association claim'})
            for f in features[p['sequence_id']]:
                fp=json.loads(f['parts_json']);dist=min(distance(a,b) for a in pp for b in fp)
                if dist>10000 or sorted(pp)==sorted(fp):continue
                k=(ck,p['record_key'],f['feature_key'])
                if k in seenpartners:continue
                seenpartners.add(k)
                partners.append({'candidate_key':ck,'sequence_id':p['sequence_id'],'rt_record_key':p['record_key'],'rt_source_branch':p['source_branch'],'neighbor_key':f['feature_key'],'neighbor_source':'OFFICIAL_CDS','neighbor_protein_accession':f['official_protein_accession'],'neighbor_product':f['product'],'neighbor_parts_json':f['parts_json'],'RT_to_neighbor_distance_bp':dist,'functional_partner_status':'NEIGHBOR_ONLY_UNPROVEN'})
            # Predicted neighbors are a distinct sensitivity layer, not official annotations.
            for f in orfby[p['sequence_id']]:
                fp=[[int(f['start0']),int(f['end0'])]];dist=min(distance(a,fp[0]) for a in pp)
                if dist>10000 or sorted(pp)==sorted(fp):continue
                k=(ck,p['record_key'],f['record_key'])
                if k in seenpartners:continue
                seenpartners.add(k)
                partners.append({'candidate_key':ck,'sequence_id':p['sequence_id'],'rt_record_key':p['record_key'],'rt_source_branch':p['source_branch'],'neighbor_key':f['record_key'],'neighbor_source':f['source_branch'],'neighbor_protein_accession':'','neighbor_product':'Prodigal predicted CDS','neighbor_parts_json':canonical_json(fp),'RT_to_neighbor_distance_bp':dist,'functional_partner_status':'NEIGHBOR_ONLY_UNPROVEN'})
    assert len(evidence)==len(candidates)==383
    for a,b in zip(candidates,evidence):
        if any(b[k]!=v for k,v in a.items()):raise AssertionError('Original candidate data changed')
    write_tsv(root/'outputs/ALL_383_RT_ORF_REPEAT_EVIDENCE.tsv',evidence,list(evidence[0]))
    write_tsv(root/'outputs/RT_CORE_GA_CANDIDATE_VIEW.tsv',[r for r in evidence if r['n_RT_core_GA_unique_AA']>0],list(evidence[0]))
    write_tsv(root/'outputs/RT_REPEAT_PAIR_CONTEXT.tsv',pairs,['candidate_key','sequence_id','rt_record_key','rt_source_branch','rt_official_accession','rt_feature_type','rt_parent_json_path','repeat_group_id','RT_to_repeat_span_distance_bp','RT_overlaps_repeat_span','repeat_high_complexity','repeat_context','repeat_predicted_CDS_overlap_bp','repeat_overlaps_original_candidate','interpretation'])
    write_tsv(root/'outputs/RT_NEIGHBOR_ORF_INVENTORY.tsv',partners,['candidate_key','sequence_id','rt_record_key','rt_source_branch','neighbor_key','neighbor_source','neighbor_protein_accession','neighbor_product','neighbor_parts_json','RT_to_neighbor_distance_bp','functional_partner_status'])
    gate={'stage':'Stage04B-R1D2','gate':'PASS_ALL383_RT_ORF_EXACT_REPEAT_EVIDENCE_COMPUTED_WITH_SCOPE_LIMITS','gate_pass':True,'n_candidates_retained':len(evidence),'n_accessions':174,'n_official_protein_records_scanned':sum(r['source_branch']=='OFFICIAL' for r in fan),'n_predicted_neighborhood_proteins_scanned':sum(r['source_branch']!='OFFICIAL' for r in fan),'n_unique_scan_targets':len({r['target_id'] for r in fan}),'n_whole_genomes_ORF_predicted':174,'n_all174_predicted_genes':len(all_orfs),'n_genome_only_rescue_predicted_genes':sum(r['sequence_id']=='NC_138483.1' for r in all_orfs),'n_candidates_with_official_RT_core_GA':sum(r['n_official_RT_core_GA']>0 for r in evidence),'n_candidates_with_any_RT_core_GA':sum(r['n_RT_core_GA_unique_AA']>0 for r in evidence),'evidence_label_counts':dict(tier_counts),'n_exact_core_repeat_groups':len(repeats),'n_candidates_with_exact_core_repeat_groups':sum(r['n_exact_core_repeat_groups']>0 for r in evidence),'profile_panel_is_exhaustive':False,'imperfect_repeat_search_exhaustive':False,'ncRNA_expression_proven':False,'ART_identities_confirmed_in_this_stage':0,'ART_identity_classification_complete':False,'novelty_audit_complete':False,'negative_calls_authorized':False,'official_missing_annotation_accessions':['NC_138483.1'],'original_official_proteins_retranslated':False,'mathematical_scores_recomputed':False,'candidate_filtering_performed':False,'rank_recomputed':False,'UG27_touched':False}
    write_json(root/'outputs/SUMMARY_STAGE04B_R1D2.json',gate)
    return gate
