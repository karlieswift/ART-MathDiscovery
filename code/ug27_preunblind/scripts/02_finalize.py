#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json, platform, sys
from pathlib import Path
import pandas as pd

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def overlap(a0,a1,b0,b1): return int(a0)<int(b1) and int(a1)>int(b0)

ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--r2c',required=True); a=ap.parse_args()
root=Path(a.root).resolve(); r2c=Path(a.r2c).resolve(); out=root/'outputs'
pre=json.loads((out/'STAGE04C2_R2D_PREFLIGHT_GATE.json').read_text())
man=json.loads((out/'CLEAN_TWO_ARM_BOARD_MANIFEST.json').read_text())
if not (pre.get('gate_pass') and man.get('gate_pass')): raise SystemExit('HOLD prior gate failed')

clean=pd.read_csv(out/'CLEAN_TWO_ARM_DUAL_ARM_VIEW.tsv',sep='\t')
legacy=pd.read_csv(r2c/'outputs/TWO_ARM_REGION_VIEW.tsv',sep='\t')

# Score-only comparison. No truth/cassette coordinates are read.
rows=[]
for _,x in clean.iterrows():
    L=legacy.loc[legacy.scoring_id.astype(str)==str(x.scoring_id)]
    ex=L.loc[(L.start0.astype(int)==int(x.start0))&(L.end0.astype(int)==int(x.end0))]
    ov=L.loc[(L.start0.astype(int)<int(x.end0))&(L.end0.astype(int)>int(x.start0))]
    rows.append({
      'assembly_id':x.assembly_id,'scoring_id':x.scoring_id,'clean_region_rank':int(x.region_rank),
      'clean_region_id':x.region_id,'clean_start0':int(x.start0),'clean_end0':int(x.end0),
      'legacy_exact_region_n':int(len(ex)),'legacy_overlap_region_n':int(len(ov)),
    })
cmp=pd.DataFrame(rows)
cmpath=out/'CLEAN_VS_LEGACY_DUAL_ARM_SCORE_ONLY_COMPARISON.tsv'; cmp.to_csv(cmpath,sep='\t',index=False)

# Output manifest before any truth is opened.
files=[p for p in out.rglob('*') if p.is_file() and p.name!='STAGE04C2_R2D_SHA256.tsv']
with (out/'STAGE04C2_R2D_SHA256.tsv').open('w') as f:
    f.write('sha256\tpath\n')
    for p in sorted(files): f.write(f'{sha(p)}\t{p.relative_to(root)}\n')

clean_counts=clean.groupby('scoring_id').size().to_dict() if len(clean) else {}
legacy_counts=legacy.groupby('scoring_id').size().to_dict() if len(legacy) else {}
final={
 'stage':'Stage04C2-R2D',
 'gate':'PASS_STAGE04C2_R2D_PREUNBLIND_CLEAN_ACNV_RECURRENCE_BOARD_FROZEN',
 'gate_pass':True,
 'n_scoring_contigs':15,
 'n_clean_union_regions':int(man['n_all_clean_regions']),
 'n_clean_dual_arm_regions':int(man['n_clean_dual_arm_regions']),
 'n_legacy_all_family_parent_dual_arm_regions':int(len(legacy)),
 'clean_dual_arm_counts_by_contig':{k:int(v) for k,v in sorted(clean_counts.items())},
 'legacy_dual_arm_counts_by_contig':{k:int(v) for k,v in sorted(legacy_counts.items())},
 'clean_exactly_equal_to_legacy_region_n':int((cmp.legacy_exact_region_n>0).sum()) if len(cmp) else 0,
 'clean_overlaps_legacy_region_n':int((cmp.legacy_overlap_region_n>0).sum()) if len(cmp) else 0,
 'families':['ACNV_FAMILY','GENERIC_RECURRENCE'],
 'proposal_q_threshold':0.005,
 'merge_rule':'strict genomic overlap; gap=0',
 'rank_rule':'support-family count desc -> sum two-family surprisal desc -> best single-family q asc -> start coordinate asc',
 'truth_loaded':False,
 'truth_paths_read':False,
 'cassette_overlap_computed':False,
 'retuning_authorized':False,
 'score_recomputed':False,
 'q_recomputed':False,
 'sequence_features_recomputed':False,
 'primary_r2c_outputs_modified':False,
 'science_change':False,
 'clean_all_regions_sha256':sha(out/'CLEAN_TWO_ARM_ALL_REGIONS.tsv'),
 'clean_dual_arm_view_sha256':sha(out/'CLEAN_TWO_ARM_DUAL_ARM_VIEW.tsv'),
 'clean_vs_legacy_sha256':sha(cmpath),
 'output_manifest_sha256':sha(out/'STAGE04C2_R2D_SHA256.tsv'),
 'next_requirement':'Archive/hash this R2D package. Only then perform a separate post-score-freeze UG27 cassette/locus unblinding without changing detector rules.',
 'runtime':{'python':sys.version,'platform':platform.platform()},
}
(out/'STAGE04C2_R2D_GATE.json').write_text(json.dumps(final,indent=2,sort_keys=True)+'\n')
# Append final gate itself to a second immutable closeout manifest.
with (out/'STAGE04C2_R2D_FINAL_SHA256.tsv').open('w') as f:
    f.write('sha256\tpath\n')
    for p in sorted([x for x in out.rglob('*') if x.is_file() and x.name!='STAGE04C2_R2D_FINAL_SHA256.tsv']):
        f.write(f'{sha(p)}\t{p.relative_to(root)}\n')
print(json.dumps(final,indent=2,sort_keys=True))
