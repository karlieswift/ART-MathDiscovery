#!/usr/bin/env python3
from collections import defaultdict

SHORTLIST_N=32
QUOTA=8


def repkey(c):
    # Byte-for-byte logical copy of historical F3B3A lines 168-175.
    return (
        0 if c['motif_c_like_in_best_domain'] else 1,
        -int(c['formal_model_n']),
        -int(c['aa_length']),
        float(c['best_domain_i_evalue']),
        c['candidate_id'],
    )


def select_families(famrows, shortlist_n=SHORTLIST_N, quota=QUOTA):
    # Exact historical F3B3A four-axis priority logic and deterministic backfill.
    if shortlist_n != 32 or quota != 8:
        raise RuntimeError('Historical F3B3A rule requires shortlist_n=32 and quota=8')
    fam_by_id={x['family_rep_id']:x for x in famrows}
    selected=[]
    reason=defaultdict(list)

    def add(items,label):
        for x in items[:quota]:
            fr=x['family_rep_id']
            reason[fr].append(label)
            if fr not in selected:
                selected.append(fr)

    add(sorted(famrows,key=lambda x:(
        -int(x['accession_n']),
        -int(x['candidate_instance_n']),
        -int(x['max_aa_length']),
        x['family_rep_id'])), 'REPLICATION_TOP8')

    add(sorted(famrows,key=lambda x:(
        -int(x['max_aa_length']),
        -int(x['accession_n']),
        -int(x['max_formal_model_n']),
        x['family_rep_id'])), 'LENGTH_TOP8')

    add(sorted(famrows,key=lambda x:(
        -int(x['max_formal_model_n']),
        -int(x['max_qualifying_model_n']),
        -int(x['accession_n']),
        -int(x['max_aa_length']),
        x['family_rep_id'])), 'MODEL_COVERAGE_TOP8')

    motif_sorted=[x for x in famrows if bool(x['motif_c_like_family'])]
    motif_sorted=sorted(motif_sorted,key=lambda x:(
        -int(x['motif_c_like_candidate_n']),
        float(x['priority_candidate_best_domain_i_evalue']),
        -int(x['priority_candidate_aa_length']),
        -int(x['accession_n']),
        x['family_rep_id']))
    add(motif_sorted, 'MOTIF_TOPOLOGY_TOP8')

    backfill=sorted(famrows,key=lambda x:(
        0 if bool(x['motif_c_like_family']) else 1,
        -int(x['accession_n']),
        -int(x['max_formal_model_n']),
        -int(x['max_aa_length']),
        -int(x['candidate_instance_n']),
        float(x['priority_candidate_best_domain_i_evalue']),
        x['family_rep_id']))
    for x in backfill:
        if len(selected) >= shortlist_n:
            break
        fr=x['family_rep_id']
        if fr not in selected:
            selected.append(fr)
            reason[fr].append('DETERMINISTIC_BACKFILL')

    if len(selected) != shortlist_n:
        raise RuntimeError(f'Historical shortlist family n={len(selected)} expected={shortlist_n}')

    rank={fr:i+1 for i,fr in enumerate(selected)}
    return selected, reason, rank, fam_by_id


def qualifying_profile_pass(definition, full_e, ie):
    if definition == 'BOTH_1E2_INCLUSIVE':
        return full_e <= 1e-2 and ie <= 1e-2
    if definition == 'IE_1E2_INCLUSIVE':
        return ie <= 1e-2
    if definition == 'FULL_1E2_INCLUSIVE':
        return full_e <= 1e-2
    if definition == 'BOTH_1E2_NONFORMAL_ONLY':
        return full_e <= 1e-2 and ie <= 1e-2 and not (full_e <= 1e-4 and ie <= 1e-4)
    if definition == 'IE_1E2_NONFORMAL_ONLY':
        return ie <= 1e-2 and not (full_e <= 1e-4 and ie <= 1e-4)
    if definition == 'FULL_1E2_NONFORMAL_ONLY':
        return full_e <= 1e-2 and not (full_e <= 1e-4 and ie <= 1e-4)
    raise RuntimeError('Unknown qualifying definition: '+definition)
