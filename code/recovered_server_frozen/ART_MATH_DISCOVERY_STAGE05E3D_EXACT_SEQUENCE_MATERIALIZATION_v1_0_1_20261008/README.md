# Stage05E3D exact sequence materialization v1.0.1

This stage closes the sequence-provenance gap exposed by E3C (`family_linked_*_sequence_n = 0`). It joins the frozen E2A candidate/witness tables to the E3A official-context/protein outputs and emits per-family exact repeat/protein FASTA files.

It does **not** change family membership, raw-DNA scores, ranks, thresholds, or E3C tiers.

## Run

```bash
./tests/run_selftest.sh
./run_stage05e3d.sh \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05E3C_FAMILY_CLOSURE_KNOWN_SYSTEM_AUDIT_v1_0_1_20261008/ART_STAGE05E3C_FAMILY_CLOSURE_RESULTS_20261008_103411.tar.gz \
  /path/to/ART_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_RESULTS_R1_20261007_172658.tar.gz \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cat outputs/STAGE05E3D_SUMMARY.json
cat outputs/STAGE05E3D_FAMILY_SEQUENCE_SUMMARY.tsv
./collect_stage05e3d.sh
```

PASS requires exact repeat + protein sequence recovery for all three frozen E3C1 families: ARCH250_0010, ARCH250_0013, ARCH250_0033. If schema drift prevents a join, the stage HOLDs and writes `STAGE05E3D_SCHEMA_DIAGNOSTICS.tsv` rather than guessing.


## v1.0.1 path hotfix
If an input path is stale because the collector wrote the result archive inside its stage directory, the runner performs an exact-basename recursive lookup under `project_root`. It accepts the recovery only when exactly one match exists; zero or multiple matches fail closed. No scientific logic changed.
