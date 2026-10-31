ART MathDiscovery Stage05E3C v1.0.1 (2026-10-08 packaging hotfix)
Family closure + known-system evidence audit + ART95 catalogue-readiness audit

HOTFIX SCOPE
- Scientific logic, E2A/E3B membership, raw-DNA scores/ranks, and E3C gates are unchanged from v1.0.0.
- Fixes the README/example E3B result path error seen on 2026-10-08.
- run_stage05e3c.sh now resolves a missing input by exact basename below PROJECT_ROOT only when there is exactly one match; otherwise it fails closed.
- Selftest no longer prints the synthetic 2-family fixture summary, avoiding confusion with a real E3C run.

PRIMARY RUN
  conda activate art_novelty
  P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
  E3C="$P/ART_MATH_DISCOVERY_STAGE05E3C_FAMILY_CLOSURE_KNOWN_SYSTEM_AUDIT_v1_0_1_20261008"
  E3B="$P/ART_MATH_DISCOVERY_STAGE05E3B_FAMILY_DEREPLICATION_v1_0_0_20261007/ART_STAGE05E3B_FAMILY_DEREPLICATION_RESULTS_20261007_182912.tar.gz"
  E3A="$P/ART_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_RESULTS_R1_20261007_172658.tar.gz"
  cd "$E3C"

  ./tests/run_selftest.sh

  ./run_stage05e3c.sh "$E3B" "$E3A" "$P"

  cat outputs/STAGE05E3C_SUMMARY.json
  column -t -s $'\t' outputs/STAGE05E3C_REVIEW_QUEUE.tsv | head -30
  cat outputs/STAGE05E3C_ART95_AUDIT_STATUS.json
  ./collect_stage05e3c.sh

NOTES
- The JSON with n_families=2 seen in the v1.0.0 selftest was generated from a deliberate two-family synthetic fixture; it was not a real-data E3C result.
- A real E3C run must observe exactly 12 frozen E3B families.
- Positive known-system evidence can trigger conflict review. Negative keyword/HMM/structure evidence is not novelty proof.
- ART95 exact 95/95 validation remains gated on an exact 95-member sequence/locus catalogue.
