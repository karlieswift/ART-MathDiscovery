Stage05F3V2 Stage04A - 60G Saturation Census

Purpose:
1) locate the frozen 6.12-Gb exact-AA unique FASTA (10,329 seqs)
2) use the frozen 54.378-Gb delta exact-AA unique FASTA (85,503 seqs)
3) globally exact-deduplicate old+delta sequences
4) re-cluster the combined exact-AA universe with the SAME frozen family rule
5) quantify exact overlap, old-only/delta-only/mixed operational families, and family-yield saturation

This is analysis-only. It does NOT rewrite old/delta family IDs, does NOT retune thresholds, and does NOT reselect the five previous systems.

Run:
  conda activate drhu
  cd /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
  tar -xzf ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04A_60G_SATURATION_CENSUS_v1_0_0_20261016.tar.gz
  cd ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04A_60G_SATURATION_CENSUS_v1_0_0_20261016
  bash run_stage04a_60g_saturation.sh /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1 64 2>&1 | tee run_stage04a.log
