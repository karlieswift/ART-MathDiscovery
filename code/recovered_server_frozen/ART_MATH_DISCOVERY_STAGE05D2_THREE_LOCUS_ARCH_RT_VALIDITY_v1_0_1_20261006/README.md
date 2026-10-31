# ART MathDiscovery Stage05D2 — Three-locus architecture + RT-validity audit

Purpose: focused post-freeze biological interpretation of three DNA-first novelty leads from Stage05D1.

Frozen loci:
- NC_031007.1::AR0009
- NC_048173.1::AR0009
- NC_062581.1::AR0001

Principle: discovery remains DNA-only (ACNV + generic recurrence). Protein evidence is post-discovery interpretation only and MUST NOT modify mathematical scores/ranks/candidate membership.

Primary questions:
1. What is the exact repeat–candidate–RT-like-ORF geometry for each locus?
2. Does the myRT broad-HMM hit have sufficient domain coverage and a recognizable RT catalytic YxDD-like motif?
3. Does current official NCBI annotation support the predicted ORF/strand/product, or is the ORF only Prodigal-predicted?
4. Do NC_031007.1 and NC_048173.1 represent duplicated/near-identical instances of one architecture family?

Run:

    conda activate art_novelty
    ./run_selftest.sh /path/to/Stage05D1
    ./run_stage05d2.sh /path/to/Stage05D1 8
    ./collect_stage05d2.sh

Outputs are review-only. No new ART identity is asserted by this stage.

## v1.0.1 hotfix (2026-10-06)
- Fixes official-annotation retrieval robustness only; no candidate, score, rank, threshold, or D1 evidence changed.
- Validates downloaded GenBank flat files before parsing (`LOCUS`, `ORIGIN`, terminal `//`).
- Uses installed NCBI EDirect `efetch` first, then documented E-utilities `rettype=gb` fallback.
- Invalid/HTML/error payloads are recorded as `FAIL_INVALID_OR_UNAVAILABLE` instead of crashing Biopython.
- GenBank parse exceptions are recorded per accession as `OFFICIAL_PARSE_FAILED` and are non-fatal.
- Replaces deprecated `Bio.pairwise2` use with `Bio.Align.PairwiseAligner`.
