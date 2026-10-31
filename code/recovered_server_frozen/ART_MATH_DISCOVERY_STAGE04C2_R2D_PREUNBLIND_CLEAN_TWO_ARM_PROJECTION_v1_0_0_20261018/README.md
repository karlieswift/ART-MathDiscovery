# Stage04C2-R2D pre-unblinding clean two-arm projection v1.0.0

Purpose: before any UG27 locus/cassette truth overlap is evaluated, derive the paper-facing **ACNV_FAMILY + GENERIC_RECURRENCE only** region boards from the already frozen Stage04C2-R2C window-family score tables using the exact historical clean-two-arm implementation.

This stage does not recompute sequence features, q-values, or detector scores. It does not read Stage04C2B-R3 truth files. It only applies the already frozen historical clean projection to the 15 frozen UG27 scoring contigs.

Runtime inputs:
1. Stage04C2-R2C root containing the completed 15-contig frozen score outputs.
2. Historical `ART_MATH_DISCOVERY_CLEAN_TWO_ARM_REPLAY_v1_0_0_20260929` root.

Historical rule reused exactly:
- families: ACNV_FAMILY and GENERIC_RECURRENCE only
- family proposal q <= 0.005
- strict genomic overlap merge, gap=0
- retain all clean union regions; no Top-K survival cutoff
- deterministic rank: support-family count desc -> sum two-family surprisal desc -> best single-family q asc -> start coordinate asc
- exact historical `scripts/common.py::build_two_arm_board()` is imported from the historical root at runtime

Frozen Stage04C2-R2C input locks:
- UG27_SCORE_FREEZE_GATE.json SHA256 `14f29cc2a856494308ffcce9b6c0d883eb10c988592ab8ebbac0e9d39447b498`
- TWO_ARM_REGION_VIEW.tsv SHA256 `1a48f06cb2991f649502d843a0e853d011c51e45cf8ba3107a467e451d1e17f3`
- HISTORICAL_EQUIVALENCE_GATE.json SHA256 `fe30d3519b693037c6e5d6705ccf06ef7129ffbac2379be09fe717518306c2ca`
- SCORE_OUTPUT_MANIFEST.tsv SHA256 `fc5baa39e6108b87f21e8e96c4691e33d57741532fdbe01b3cd6d4955d9efc95`

Expected historical clean implementation hashes:
- scripts/common.py `a59b0e8439747487238b25249116f476baf54e55a5a6681411cb2272763cf255`
- scripts/01_build_boards.py `e0f7de502314f1854b2520b97892b56052b05b9db1f0be88c2ba15aa3493680a`
- tests/test_two_arm.py `8ff3f567e3eb6c040a105e5ca45038b9e840100365a3853430091812c7b902da`
- README.md `91ef9fef928a0969fa5393903874f5bbd3ee6bb81ac910add19246c5c8a0dee2`

Run:

```bash
./run_selftest.sh
./run_stage04c2_r2d.sh \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE04C2_R2C_FULL_STRATUM_REPLAY_v1_0_0_20261005 \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_CLEAN_TWO_ARM_REPLAY_v1_0_0_20260929
```

Primary gate:
`outputs/STAGE04C2_R2D_GATE.json`

Only after that gate passes and the output manifest is archived/hashed should UG27 cassette/locus truth be loaded for the separate post-score-freeze unblinding stage.
