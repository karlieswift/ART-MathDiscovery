Stage05F3V2 Stage04B1 - Delta-family replication census
Date: 2026-10-16

Purpose
- Read only the frozen Stage04A combined family census and frozen 54.378-Gb delta formal FASTAs.
- Reconstruct family presence across the 254 source-blind delta accessions.
- Quantify how many DELTA_ONLY and MIXED_OLD_DELTA families recur in >=2, >=3, >=5, >=10 accessions.
- If the selected manifest contains exactly one explicit *stratum* column, also report source-stratum breadth.

Scientific boundary
- No new HMM scan.
- No re-clustering.
- No threshold tuning.
- No annotation or truth used for selection.
- No structure prediction.
- No priority selection.
- Multi-accession recurrence is a replication descriptor, not proof of biological independence/function.

Run
  conda activate drhu
  P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
  cd $P/ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS_v1_0_0_20261016
  bash run_stage04b1_delta_replication.sh "$P" 2>&1 | tee run_stage04b1.log

Expected terminal gate
  PASS_STAGE05F3V2_STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS_COMPLETE
