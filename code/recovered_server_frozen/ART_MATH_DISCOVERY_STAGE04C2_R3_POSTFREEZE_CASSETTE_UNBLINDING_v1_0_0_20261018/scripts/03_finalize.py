#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,platform,sys
from pathlib import Path
import pandas as pd
from common import sha256_file

ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); a=ap.parse_args()
root=Path(a.root).resolve(); out=root/'outputs'
S=pd.read_csv(out/'UG27_CASSETTE_UNBLINDING_SUMMARY.tsv',sep='\t')
if len(S)!=5: raise SystemExit('HOLD_STAGE04C2_R3_FINALIZE: expected five systems')
classes=S['cassette_level_classification'].value_counts().to_dict()
primary=int(S['primary_cassette_dual_arm_strict_overlap'].astype(bool).sum())
D={
 'stage':'Stage04C2-R3',
 'gate':'PASS_STAGE04C2_R3_POSTFREEZE_CASSETTE_LEVEL_CLEAN_TWO_ARM_UNBLINDING_COMPLETE_NO_RETUNE',
 'gate_pass':True,
 'n_prespecified_systems':5,
 'n_evaluated_systems':5,
 'primary_endpoint':'strict 0-based half-open coordinate overlap of frozen native cassette interval with frozen clean ACNV+GENERIC_RECURRENCE dual-arm region',
 'cassette_dual_arm_strict_overlap_n':primary,
 'cassette_dual_arm_strict_overlap_fraction':primary/5.0,
 'classification_counts':{str(k):int(v) for k,v in classes.items()},
 'truth_type':'WT construct-to-native long-HSP cassette interval',
 'exact_ncrna_array_truth':False,
 'nearby_threshold_defined':False,
 'nearest_region_gap_is_descriptive_continuous_quantity':True,
 'retuning_authorized':False,
 'score_recomputed':False,
 'q_recomputed':False,
 'sequence_features_recomputed':False,
 'region_board_recomputed':False,
 'region_rank_recomputed':False,
 'science_change':False,
 'interpretation_boundary':'This is cassette-level localization against five frozen WT construct-to-native intervals. It does not establish exact ncRNA-array overlap, biochemical RT activity, ART identity, or novelty.',
 'summary_sha256':sha256_file(out/'UG27_CASSETTE_UNBLINDING_SUMMARY.tsv'),
 'direct_overlaps_sha256':sha256_file(out/'UG27_CASSETTE_DIRECT_OVERLAPS.tsv'),
 'nearest_regions_sha256':sha256_file(out/'UG27_CASSETTE_NEAREST_REGIONS.tsv'),
 'truth_lock_sha256':sha256_file(out/'UG27_R3_CASSETTE_TRUTH_LOCK.tsv'),
 'score_archive_preflight_sha256':sha256_file(out/'STAGE04C2_R3_SCORE_ARCHIVE_PREFLIGHT_GATE.json'),
 'runtime':{'python':sys.version,'platform':platform.platform()},
}
(out/'STAGE04C2_R3_GATE.json').write_text(json.dumps(D,indent=2,sort_keys=True)+'\n')
with (out/'STAGE04C2_R3_SHA256.tsv').open('w') as f:
    f.write('sha256\tpath\n')
    for p in sorted([x for x in out.rglob('*') if x.is_file() and x.name!='STAGE04C2_R3_SHA256.tsv']):
        f.write(f'{sha256_file(p)}\t{p.relative_to(root)}\n')
print(json.dumps(D,indent=2,sort_keys=True))
