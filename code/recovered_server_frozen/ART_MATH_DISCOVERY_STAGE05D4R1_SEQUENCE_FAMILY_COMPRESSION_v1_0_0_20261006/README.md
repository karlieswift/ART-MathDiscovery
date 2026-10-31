# Stage05D4-R1 sequence-family compression

Consumes the frozen Stage05D4 N1/N2 nonformal broad-RT protein universe (13,062 unique AA). It performs two fixed MMseqs2 clusterings without changing any DNA candidate, score, or rank:

- NR90: min-seq-id 0.90, coverage 0.80, cov-mode 0. Used only to collapse near-redundant queries for downstream remote-homology work.
- FAM50: min-seq-id 0.50, coverage 0.70, cov-mode 0. Descriptive family map only; never a candidate filter or identity call.

No RT, ART, novelty, or wet-lab identity is assigned in this stage.
