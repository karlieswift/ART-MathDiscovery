# Stage05D4-R5 v1.0.8 postprocessing-only hotfix

Observed after the successful v1.0.7 Foldseek run:
- Foldseek completed and wrote `raw/P1_FULL3D_AFSP.m8` with 2 rows.
- The legacy summarizer reported `n_queries_with_full3d_hit=0`.

Cause: the frozen R5 query IDs are `AA_<sha>`, while Foldseek created from Boltz structure files can emit file-derived IDs such as `AA_<sha>_model_0`. The legacy summarizer keyed hits by the raw Foldseek query string and therefore failed to join valid m8 rows back to the frozen query universe.

v1.0.8 changes only postprocessing. It reuses the existing 43 Boltz predictions and the existing v1.0.7 Foldseek m8 file; no structure prediction or Foldseek search is rerun. It normalizes each m8 query ID to exactly one frozen P1 query, writes an explicit query-ID mapping audit table, and regenerates the triage/summary/gate outputs.

Scientific settings, P1=43, AF/Swiss-Prot target, E-value, max-seqs, exact-TMscore, DNA scores/ranks, and candidate universe are unchanged.
