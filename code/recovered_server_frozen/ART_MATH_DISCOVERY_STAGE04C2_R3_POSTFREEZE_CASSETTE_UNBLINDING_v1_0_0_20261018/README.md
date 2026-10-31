# Stage04C2-R3 post-score-freeze cassette unblinding v1.0.0

Purpose: evaluate the already-frozen Stage04C2-R2D clean ACNV+GENERIC_RECURRENCE region board against the five already-frozen UG27 WT construct-to-native cassette intervals.

This stage does not recompute sequence features, q values, proposals, clean regions, region ranks, or detector parameters. It first verifies that the exact R2D result archive and live R2D outputs are byte-locked, then opens the previously frozen Stage04C2B-R3 bridge.

Primary endpoint: strict coordinate overlap between each native cassette interval and any frozen clean dual-arm region. The cassette is the native interval spanned by the long WT construct HSP. It is NOT automatically the exact ncRNA-array truth interval. Therefore all outputs are explicitly cassette-level localization results.

If no clean dual-arm region overlaps a cassette, the stage separately records direct ACNV-only and recurrence-only clean regions. It distinguishes both single arms present but region-level nonconvergent from one-arm-only direct support. Nearest-region gaps are continuous descriptive quantities; no post-hoc nearby cutoff is introduced.

Classification hierarchy:
- DUAL_ARM_DIRECT
- BOTH_SINGLE_ARMS_DIRECT_NONCONVERGENT
- ACNV_ONLY_DIRECT
- RECURRENCE_ONLY_DIRECT
- NO_DIRECT_CLEAN_REGION

No performance threshold determines gate PASS. The gate means the five prespecified cassette-level evaluations completed without changing the frozen detector.
