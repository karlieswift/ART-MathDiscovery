# ART panel StageA1: recover existing Stage05B frozen scores

The five newly frozen exact sources (PB50, LY01, S6, PALS_2, UFV_DC4) were already members of the 919-genome Stage05A/Stage05B universe. Therefore **do not rescore them merely to make an ART panel**. This package recovers the existing Stage05B frozen proposal/union evidence and binds it to the exact FASTA SHA256 values from `ART_EXACT_SOURCE_REPLAY_INPUTS_*`.

No score, feature, q threshold, lag, region rank, or detector parameter is recomputed. Array coordinates and repeat motifs remain sealed.

Run:

```bash
./tests/run_selftest.sh
./run_recover_existing_scores.sh \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1 \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_EXACT_SOURCE_REPLAY_INPUTS_20261008_111922.tar.gz

cat outputs/ART_PANEL_STAGEA1_GATE.json
cat outputs/ART_PANEL_STAGE05B_SCORE_RECOVERY_SUMMARY.tsv
./collect_results.sh
```

If the Stage05B archive is not at the project root, the runner performs an exact-basename recursive lookup for `ART_STAGE05B_NEW919_SCORE_RESULTS_20261006_081746.tar.gz` and accepts only a unique match. The expected archive SHA256 is frozen as `4fc4f1f09ed92f4e0be9817b14c0c14fe7f9eafd56cc0918cad5e6a73a054ebd`.
