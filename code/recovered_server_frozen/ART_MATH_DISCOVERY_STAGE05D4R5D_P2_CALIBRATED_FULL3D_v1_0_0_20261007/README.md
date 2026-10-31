# ART MathDiscovery Stage05D4-R5D — P2 calibrated full-3D extension

Purpose: run the already-prespecified P2 compact singleton tier (n=40) through the same MSA-assisted Boltz-2 full-coordinate workflow and the R5C-calibrated clean-RT Foldseek screen.

This is a compute-priority extension, not a new candidate filter. It does not change DNA scores/ranks, P2 membership, Boltz settings, or R5C structural thresholds.

## Inputs
- R3 results root containing `outputs/STAGE05D4R3_TIER_A_COMPACT_STRUCTURE_READY.faa`
- R5 root containing `frozen/STAGE05D4R5_P2_EXTENSION_FREEZE.tsv` and reusable `boltz_cache`
- R5C root with PASS calibration and `references/clean_rt_afsp`

## Run
1. `./run_selftest.sh`
2. `./prepare_stage05d4r5d.sh R3_RESULTS_ROOT R5_ROOT`
3. `./run_stage05d4r5d_p2.sh R5_ROOT 2`
4. `./run_clean_rt_screen.sh R5C_ROOT 32`
5. `./collect_stage05d4r5d.sh`

The Boltz step deliberately uses four independent single-GPU workers (10 queries each) and `--no_kernels`, avoiding the DDP/NCCL and cuequivariance issues already audited in R5.
