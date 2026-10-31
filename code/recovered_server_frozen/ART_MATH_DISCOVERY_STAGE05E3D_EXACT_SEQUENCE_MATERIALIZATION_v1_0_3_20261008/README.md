# ART MathDiscovery Stage05E3D v1.0.3

Corrected exact-witness sequence materialization. **This supersedes v1.0.1 for all protein-family downstream analysis.**

v1.0.1's 1,629 protein sequences were dominated by general neighborhood proteins because a permissive fallback token join was applied to the entire official-neighborhood table. The E3C family conclusions remain unchanged, but those 1,629 sequences must not be treated as family members.

v1.0.3 uses only explicit E3A witness links and requires exact concordance with the frozen E3C official-witness counts/median lengths. Priority families ARCH250_0010/0013/0033 must be sequence-complete.

Run:

```bash
./tests/run_selftest.sh
./run_stage05e3d.sh <E3C_results.tar.gz|dir> <E3A_results.tar.gz|dir> <project_root>
cat outputs/STAGE05E3D_SUMMARY.json
cat outputs/STAGE05E3D_FAMILY_SEQUENCE_SUMMARY.tsv
./collect_stage05e3d.sh
```

## v1.0.3 hotfix

Fixes an alias-map implementation bug in v1.0.2: if an unversioned FASTA alias had already become ambiguous, a third record sharing that alias could dereference `None`. v1.0.3 keeps such aliases explicitly ambiguous while preserving exact versioned IDs. No sequence-selection, family, score, rank, threshold, or concordance rule changed.
