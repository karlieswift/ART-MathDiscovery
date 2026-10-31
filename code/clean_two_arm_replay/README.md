# ART clean two-arm replay v1.0.0

Purpose: derive a paper-facing **ACNV + generic recurrence only** region board from already-frozen pre-unblinding window score tables. This package does **not** recompute ACNV/recurrence features, change q-values, alter thresholds, or use truth coordinates before the new boards are written and hashed.

Frozen rule:
- families: `ACNV_FAMILY`, `GENERIC_RECURRENCE`
- proposal: family q <= 0.005
- merge: strict genomic overlap, gap=0
- retain all union regions; no Top-K cutoff
- rank: support-family count desc -> sum two-family surprisal desc -> best single-family q asc -> coordinate asc

Inputs:
1. MarsHill Stage02B root containing `outputs/STAGE02B_WINDOW_FAMILY_SCORES.tsv.gz`.
2. Exact Stage03A results tar.gz. Expected SHA256 is `d46b0246ae4cf965b18663d27eee04ca7c917596d96a5394f44acf87abde6308`.

Recommended command on the project server:

```bash
cd /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
unzip -o ART_MATH_DISCOVERY_CLEAN_TWO_ARM_REPLAY_v1_0_0_20260929.zip
cd ART_MATH_DISCOVERY_CLEAN_TWO_ARM_REPLAY_v1_0_0_20260929
./run_selftest.sh
./run_clean_two_arm_replay.sh \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE02B_MULTI_AGENT_REGION_SUPERVISOR_v1_0_0_20260927 \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE03A_FROZEN_ENGINE_L3_VALIDATION_v1_0_0_20260927/ART_MATH_DISCOVERY_STAGE03A_RESULTS_20260928_085109.tar.gz
```

Primary outputs:
- `outputs/CLEAN_TWO_ARM_PREFLIGHT_GATE.json`
- `outputs/boards/<case>/CLEAN_TWO_ARM_UNION_REGIONS.tsv`
- `outputs/CLEAN_TWO_ARM_POSTFREEZE_UNBLINDING.tsv`
- `outputs/CLEAN_TWO_ARM_FINAL_GATE.json`
- `PACKAGE_MANIFEST.sha256`

Evidence status: a post-hoc paper-simplification replay from frozen pre-unblinding scores, not a new prospective experiment. Truth is used only after the four clean boards are built.
