# ART MathDiscovery Stage05E2 — Global specificity + topology audit

Purpose: convert the 66 Stage05E1 primary replicated architecture families into a review-order audit without changing family membership. Specificity is measured against all 999 Stage05D1 novelty-first candidates. No P1/P2 structure result is used.

Inputs:
1. Stage05D1 novelty-first package root.
2. Stage05E1 package root containing `outputs/`.
3. Stage05D4R1 sequence-family-compression results tar.gz or extracted directory.

Run:
```bash
./run_selftest.sh
./run_stage05e2.sh <D1_ROOT> <E1_ROOT> <R1_RESULTS_TAR_OR_DIR>
```

Review strata:
- E2A: exact-primary repeat-core + compact FAM50 joint family-specific signature.
- E2B: approximate/mixed repeat-core + compact FAM50 joint family-specific signature.
- E2C: partial family-specific evidence, no joint edge.
- E2D: broadly shared evidence under this frozen audit representation.

All 66 primary families and 324 candidate loci remain retained. These are review-order strata, not ART calls.
