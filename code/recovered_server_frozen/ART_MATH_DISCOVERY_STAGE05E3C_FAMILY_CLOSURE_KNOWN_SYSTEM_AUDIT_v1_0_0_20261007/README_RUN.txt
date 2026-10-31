ART MathDiscovery Stage05E3C v1.0.0
Family closure + known-system evidence audit + ART95 catalogue-readiness audit

PURPOSE
- Continue from frozen E3A/E3B only.
- Keep all E2A/E3B family membership unchanged.
- Do NOT recompute raw-DNA ACNV/recurrence scores or ranks.
- Treat positive known-system evidence as a reason for conflict review.
- Never treat a negative keyword/HMM/structure result as novelty proof.
- ART95: an exact 95/95 validation requires an exact 95-member sequence/locus catalogue. The supplied public-reference table is deliberately incomplete and cannot support a 95/95 claim.

PRIMARY RUN
  conda activate art_novelty
  P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
  E3C="$P/ART_MATH_DISCOVERY_STAGE05E3C_FAMILY_CLOSURE_KNOWN_SYSTEM_AUDIT_v1_0_0_20261007"
  cd "$E3C"

  ./tests/run_selftest.sh

  ./run_stage05e3c.sh \
    "$P/ART_STAGE05E3B_FAMILY_DEREPLICATION_RESULTS_20261007_182912.tar.gz" \
    "$P/ART_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_RESULTS_R1_20261007_172658.tar.gz" \
    "$P"

  cat outputs/STAGE05E3C_SUMMARY.json
  column -t -s $'\t' outputs/STAGE05E3C_REVIEW_QUEUE.tsv | head -30
  cat outputs/STAGE05E3C_ART95_AUDIT_STATUS.json

  ./collect_stage05e3c.sh

OPTIONAL: exact ART95 catalogue check
  export E3C_ART95_MANIFEST=/path/to/ART95_EXACT_CATALOGUE.tsv
  ./run_stage05e3c.sh <E3B> <E3A> "$P"

OPTIONAL: local HMM conflict scan
  export E3C_KNOWN_HMM_DIR=/path/to/known_system_hmms
  ./run_stage05e3c.sh <E3B> <E3A> "$P"

The optional HMM scan is POST-SCORE interpretation only. A hit can flag a family for known-system review; a no-hit cannot establish novelty.
