# Stage05D4-R5 v1.0.5 Foldseek query-DB hotfix

Purpose: downstream-only compatibility repair after 43/43 Boltz-2 full-coordinate predictions are already complete.

Observed state before this hotfix:
- canonical `boltz_out/predictions/` contains 43 real query directories;
- 43 regular `*_model_0.cif` files and 43 regular confidence JSON files are present;
- previous full3D Foldseek wrapper flattened CIFs as symbolic links, then passed that directory directly to `foldseek easy-search`;
- Foldseek 10.941cd33 returned after parameter printing without producing the expected downstream summary.

v1.0.5 removes that filesystem ambiguity by:
1. materializing the 43 query CIFs into a fresh flat directory as **regular files** (hard-link preferred, copy fallback);
2. building an explicit reusable Foldseek query structure DB with `foldseek createdb`;
3. verifying the created query DB contains exactly 43 entries via `convert2fasta`;
4. searching the explicit query DB against the already frozen AlphaFold/Swiss-Prot target DB;
5. running the existing frozen Stage05D4-R5 summarizer unchanged.

Scientific settings are unchanged. This hotfix does not rerun Boltz, alter P1, change sequence/MSA/model/seed, recompute DNA scores/ranks, change the Foldseek target database, or change the E-value / max-seqs / exact-TMscore search settings.

Also note: earlier v1.0.4 `.sha256` sidecar accidentally embedded an absolute `/mnt/data/...` source path. That was a packaging-only sidecar bug. This v1.0.5 sidecar uses a relative basename and can be checked from the project root with `sha256sum -c`.
