"""Exact-core direct-repeat evidence. Does not infer RNA expression or ART identity."""
from __future__ import annotations
from collections import defaultdict,Counter
from math import log2
import json
from common import *

def entropy(s):
    return -sum((n/len(s))*log2(n/len(s)) for n in Counter(s).values()) if s else 0.0

def find_exact_core_groups(seq,policy):
    k=policy['seed_nt'];minc=policy['min_copies'];maxcore=policy['max_core_nt'];maxstep=policy['max_adjacent_start_distance_bp']
    occur=defaultdict(list)
    for i in range(len(seq)-k+1):
        motif=seq[i:i+k]
        if not (set(motif)-set('ACGT')): occur[motif].append(i)
    candidates=[]
    for motif,ps in sorted(occur.items()):
        if len(ps)<minc:continue
        # Connected components on adjacent start distances. An exact seed is
        # literal on the forward sequence; the same direct repeats also exist
        # on the reverse complement, so no second count is needed.
        clusters=[];cluster=[]
        for p in ps:
            if cluster and p-cluster[-1]>maxstep:clusters.append(cluster);cluster=[]
            cluster.append(p)
        if cluster:clusters.append(cluster)
        for ps in clusters:
            chosen=[]
            for p in ps:
                if not chosen or p-chosen[-1]>=k:chosen.append(p)
            if len(chosen)<minc:continue
            gaps=[b-a for a,b in zip(chosen,chosen[1:])]
            # Keep the greedy non-overlapping path; never silently claim this
            # is an exhaustive approximate-repeat or periodicity algorithm.
            if max(gaps)>maxstep:continue
            cap=min(maxcore,min(gaps));left=0;right=k
            while right-left<cap and min(chosen)+left>0:
                chars={seq[p+left-1] for p in chosen}
                if len(chars)!=1:break
                left-=1
            while right-left<cap and max(chosen)+right<len(seq):
                chars={seq[p+right] for p in chosen}
                if len(chars)!=1:break
                right+=1
            starts=tuple(p+left for p in chosen);core=seq[starts[0]:chosen[0]+right];L=len(core)
            if any(seq[p:p+L]!=core for p in starts):raise AssertionError('Nonexact repeat emitted')
            candidates.append({'starts':starts,'repeat_len':L,'core':core,'seed_count':1})
    # Merge only representations whose copies are contained one-for-one.
    # Do not merge differing copy counts by a subjective best-array score.
    kept=[];seen={}
    for c in sorted(candidates,key=lambda c:(-c['repeat_len'],-len(c['starts']),c['starts'],c['core'])):
        key=(c['starts'],c['repeat_len'])
        if key in seen:seen[key]['seed_count']+=1;continue
        match=None
        for old in kept:
            if len(old['starts'])!=len(c['starts']):continue
            if all(a<=b and b+c['repeat_len']<=a+old['repeat_len'] for a,b in zip(old['starts'],c['starts'])):match=old;break
        if match:match['seed_count']+=1
        else:kept.append(c);seen[key]=c
    return sorted(kept,key=lambda c:(c['starts'][0],c['starts'][-1]+c['repeat_len'],c['core']))

REPEAT_FIELDS=['candidate_key','sequence_id','repeat_group_id','start0','end0','repeat_core_nt','core_length_nt','n_copies','copy_intervals_json','spacer_lengths_json','spacing_cv','core_entropy_bits','high_complexity_evidence','n_seed_representations','overlaps_frozen_candidate','distance_to_candidate_bp','official_annotation_available','official_cds_overlap_bp_span','official_cds_overlap_bp_copies','n_copies_overlapping_official_cds','annotation_relative_repeat_context','predicted_orf_assessed','predicted_cds_overlap_bp_span','predicted_cds_overlap_bp_copies','rna_expression_proven','ART_identity_call']
WINDOW_FIELDS=['candidate_key','sequence_id','candidate_start0','candidate_end0','recurrence_start0','recurrence_end0','official_annotation_available','candidate_official_cds_overlap_bp','recurrence_official_cds_overlap_bp','recurrence_official_cds_overlap_fraction','predicted_orf_assessed','recurrence_predicted_cds_overlap_bp','n_exact_core_repeat_groups','n_high_complexity_repeat_groups','n_high_complexity_full_span_outside_official_cds','n_high_complexity_full_span_outside_official_and_predicted_cds','interpretation']

def build_repeat_evidence(root,predicted=None):
    root=Path(root); pol=read_json(root/'config/METHOD_POLICY.json'); cand=read_tsv(root/'inputs/candidates.tsv')
    sources={r['accession']:r for r in read_tsv(root/'inputs/genome_sources.tsv')}
    genomes={acc:next(fasta(root/'inputs/genomes'/f'{acc}.fna'))[2] for acc in sources}
    official=defaultdict(list)
    for f in read_tsv(root/'inputs/all_official_features.tsv'):
        if f['feature_type']=='CDS':official[f['sequence_id']].extend(json.loads(f['parts_json']))
    masks={acc:merge_intervals(official[acc]) for acc in sources}
    pmm=defaultdict(list)
    if predicted is not None:
        for r in predicted:pmm[r['sequence_id']].append((int(r['start0']),int(r['end0'])))
    pmasks={acc:merge_intervals(pmm[acc]) for acc in sources}
    reps=[];windows=[]
    for idx,c in enumerate(cand,1):
        acc=c['sequence_id'];key=acc+'::'+c['region_id'];L=len(genomes[acc]);c0,c1=int(c['start0']),int(c['end0']);n0=max(0,c0-pol['neighborhood_bp']);n1=min(L,c1+pol['neighborhood_bp'])
        annot=sources[acc]['annotation_status']=='OFFICIAL_CDS_AVAILABLE'; mask=masks[acc]; pm=pmasks[acc]
        found=find_exact_core_groups(genomes[acc][n0:n1],pol['repeat_survey']);local=[]
        for i,r in enumerate(found,1):
            parts=[[n0+p,n0+p+r['repeat_len']] for p in r['starts']];span=(parts[0][0],parts[-1][1]); steps=[b[0]-a[0] for a,b in zip(parts,parts[1:])];m=sum(steps)/len(steps);cv=(sum((s-m)**2 for s in steps)/len(steps))**0.5/m
            cover=covered_bp(span,mask);cc=[covered_bp(x,mask) for x in parts]
            context='ANNOTATION_ABSENT' if not annot else ('FULL_SPAN_OUTSIDE_ANNOTATED_CDS' if cover==0 else ('COPIES_OUTSIDE_BUT_SPACERS_OVERLAP_CDS' if sum(cc)==0 else 'CODING_OR_MIXED'))
            ent=entropy(r['core']);high=ent>=pol['repeat_survey']['complexity_high_entropy_bits_min']
            out={'candidate_key':key,'sequence_id':acc,'repeat_group_id':key+f'::DREP{i:04d}','start0':span[0],'end0':span[1],'repeat_core_nt':r['core'],'core_length_nt':r['repeat_len'],'n_copies':len(parts),'copy_intervals_json':canonical_json(parts),'spacer_lengths_json':canonical_json([b[0]-a[1] for a,b in zip(parts,parts[1:])]),'spacing_cv':cv,'core_entropy_bits':ent,'high_complexity_evidence':int(high),'n_seed_representations':r['seed_count'],'overlaps_frozen_candidate':int(overlaps(span,(c0,c1))),'distance_to_candidate_bp':distance(span,(c0,c1)),'official_annotation_available':int(annot),'official_cds_overlap_bp_span':cover if annot else '', 'official_cds_overlap_bp_copies':sum(cc) if annot else '', 'n_copies_overlapping_official_cds':sum(x>0 for x in cc) if annot else '', 'annotation_relative_repeat_context':context,'predicted_orf_assessed':int(predicted is not None),'predicted_cds_overlap_bp_span':covered_bp(span,pm) if predicted is not None else '', 'predicted_cds_overlap_bp_copies':sum(covered_bp(x,pm) for x in parts) if predicted is not None else '', 'rna_expression_proven':False,'ART_identity_call':'NOT_ASSESSED'}
            local.append(out);reps.append(out)
        rw=(int(c['recurrence_best_window_start0']),int(c['recurrence_best_window_end0']))
        windows.append({'candidate_key':key,'sequence_id':acc,'candidate_start0':c0,'candidate_end0':c1,'recurrence_start0':rw[0],'recurrence_end0':rw[1],'official_annotation_available':int(annot),'candidate_official_cds_overlap_bp':covered_bp((c0,c1),mask) if annot else '', 'recurrence_official_cds_overlap_bp':covered_bp(rw,mask) if annot else '', 'recurrence_official_cds_overlap_fraction':covered_bp(rw,mask)/(rw[1]-rw[0]) if annot else '', 'predicted_orf_assessed':int(predicted is not None),'recurrence_predicted_cds_overlap_bp':covered_bp(rw,pm) if predicted is not None else '', 'n_exact_core_repeat_groups':len(local),'n_high_complexity_repeat_groups':sum(r['high_complexity_evidence'] for r in local),'n_high_complexity_full_span_outside_official_cds':sum(r['high_complexity_evidence'] and r['annotation_relative_repeat_context']=='FULL_SPAN_OUTSIDE_ANNOTATED_CDS' for r in local) if annot else '', 'n_high_complexity_full_span_outside_official_and_predicted_cds':sum(r['high_complexity_evidence'] and r['annotation_relative_repeat_context']=='FULL_SPAN_OUTSIDE_ANNOTATED_CDS' and r['predicted_cds_overlap_bp_span']==0 for r in local) if predicted is not None and annot else '', 'interpretation':'Descriptive exact direct-repeat context, not RNA-expression, noncoding, or ART proof'})
        if idx%25==0:print(f'[REPEATS] {idx}/{len(cand)} groups={len(reps)}',flush=True)
    return reps,windows
