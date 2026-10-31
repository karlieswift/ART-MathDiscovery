#!/usr/bin/env python3
import csv, json, sys
from pathlib import Path
b7=Path(sys.argv[1]).resolve()
p=b7/'outputs/MIXED23_FAMILYWIDE_HISTORICAL_CROSSWALK.tsv'
rows=list(csv.DictReader(open(p),delimiter='\t'))
fan=[]
for r in rows:
    h=[x for x in r.get('historical5_overlap','').split(';') if x]
    if len(h)>=2:
        fan.append((r['delta_operational_family_id'],len(h),h,int(r.get('direct_pair_hit_n','0') or 0),int(r.get('mapped_old_family_n','0') or 0)))
all5=sorted({x for r in rows for x in r.get('historical5_overlap','').split(';') if x})
out={
 'stage':'Stage05F3V2-Stage04B7R1',
 'historical5_with_any_pairwise_connection_n':len(all5),
 'mixed23_with_multi_historical5_fanout_n':len(fan),
 'specific_system_replication_claim_authorized':len(fan)==0,
 'authorized_interpretation':('ONE_TO_ONE_SPECIFIC_REPLICATION_NOT_DISQUALIFIED_BY_THIS_GUARD' if len(fan)==0 else 'PAIRWISE_SEQUENCE_SPACE_COMPATIBILITY_ONLY_PENDING_FROZEN_COMBINED_CLUSTER_CROSSWALK'),
 'retuning_performed':False,
 'science_threshold_changed':False,
}
print(json.dumps(out,indent=2))
if fan:
    print('\nMULTI-HISTORICAL5 FANOUT:')
    for fam,n,hs,npair,nold in fan:
        print(f'{fam}\thistorical5_n={n}\tpair_hits={npair}\tmapped_old_family_n={nold}\t'+';'.join(hs))
    print('\nHOLD_SPECIFIC_HISTORICAL5_SYSTEM_REPLICATION_INTERPRETATION_MANY_TO_MANY')
else:
    print('PASS_B7_REPLICATION_SEMANTIC_GUARD_NO_MULTI_H5_FANOUT')
