# Stage05D4-R5 v1.0.6 downstream-only hotfix

Diagnosis:
- v1.0.5 successfully built a 43-entry Foldseek query database.
- But v1.0.5 then passed that query DB prefix to `foldseek easy-search`.
- Foldseek's easy-search workflow accepts raw query structure file(s)/folder(s); a prebuilt target DB is supported. For prebuilt queryDB + targetDB the lower-level workflow is `foldseek search` followed by conversion.
- The user's manual rerun also occurred from the project root, so relative `work/...` paths referred to the wrong directory.

Fix:
- `cd` into R5 explicitly.
- Reuse the already verified 43 regular CIF files in `work/p1_cif_flat_regular`.
- Run `foldseek easy-search <regular-CIF-dir> <AF/Swiss-Prot DB> ...` with the unchanged frozen search parameters.
- No Boltz prediction rerun, no DNA rescoring/reranking, no candidate filtering, no scientific threshold/model changes.
