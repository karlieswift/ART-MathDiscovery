# ART MathDiscovery Stage05E3D v1.0.5 — predicted-witness hash rescue

This release changes the sequence-retrieval target, not the frozen scientific universe.

## Why v1.0.4 failed

E3A's `official-supported` classification is an annotation/coordinate support layer. It does **not** guarantee that every `witness_official_protein_id` is present in the exported `STAGE05E3A_NEIGHBORHOOD_OFFICIAL_PROTEINS.faa`. The real server error for `KY555142.1::AR0003 / ARB13856.1` proved that assumption false (`pair=absent; exact_id=absent; unversioned=absent`).

More importantly, the E1/E2 family-defining protein identity is the **predicted witness ORF/FAM50 member**, frozen in `STAGE05E3A_E2A_WITNESS_PAIRS.tsv` by `record_key` + `aa_sha256`. The official protein is a support annotation, not the sequence that defined family membership.

v1.0.5 therefore:

1. keeps the E3C 12-family universe and E3A official-support labels unchanged;
2. selects the exact E2A predicted witness for each E3A official-supported candidate;
3. recovers its amino-acid sequence by exact SHA256 from existing Stage05 protein FASTAs under the project root;
4. materializes exact E2A repeat cores;
5. never substitutes a neighboring protein and never requires the official protein FASTA to contain the annotation ID.

No DNA score/rank, family membership, threshold, ART identity, novelty, function, or wet-lab claim is changed.

## Run

```bash
./tests/run_selftest.sh
./run_stage05e3d.sh <E3C_results.tar.gz|dir> <E3A_results.tar.gz|dir> <project_root>
cat outputs/STAGE05E3D_SUMMARY.json
column -t -s $'\t' outputs/STAGE05E3D_FAMILY_SEQUENCE_SUMMARY.tsv
cat outputs/STAGE05E3D_WITNESS_SELECTION_DIAGNOSTICS.tsv | head -30
cat outputs/STAGE05E3D_MISSING_PREDICTED_WITNESS_SEQUENCES.tsv
./collect_stage05e3d.sh
```

Expected full PASS requires all E3C official-supported witnesses to resolve to exact E2A hashes and all those hashes to be found in frozen project protein FASTAs. If a candidate has more than one possible E2A hash and no exact record/hash link in E3A context, the stage holds instead of guessing.
