# Stage05E3F — priority architecture expansion + RT-context closure v1.0.1

Purpose: read the already-PASS Stage05E3E result, separate protein-family expansion from exact repeat-linked architecture expansion for ARCH250_0010 / 0013 / 0033, and close Stage05C1 RT context without rescoring anything.

Hotfix relative to v1.0.0: the original package required `RT_LOCI_ALL919.tsv` to still be present as a live file or inside the Stage05C1 result tar. On the production project tree neither resource was available. v1.0.1 therefore keeps the exact-coordinate route when available, but adds a fail-closed frozen accession-level route using the Stage05C1/C2A/C2B ten-RT-core-genome ledger. Accession-level closure is permitted ONLY if none of the E3E homolog genomes intersects the ten frozen RT-core genomes. If any accession intersects, exact RT coordinates are required and the run stops.

For the received E3E result the local homolog accessions are disjoint from the frozen ten-RT-genome set, so a zero same-genome RT-context conclusion does not require RT coordinates. This does not mean the proteins are globally non-RT-related; it only means none occurs in a genome carrying one of the ten Stage05C1 RT-core GA loci.

Run:
```bash
conda activate art_rt_evidence
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P/ART_MATH_DISCOVERY_STAGE05E3F_PRIORITY_ARCHITECTURE_EXPANSION_RT_CONTEXT_CLOSURE_v1_0_1_20261008"
./tests/run_selftest.sh
./run_stage05e3f.sh \
  "$P/ART_STAGE05E3E_PRIORITY_FAMILY_BIOLOGICAL_CLOSURE_RESULTS_20261008_141923.tar.gz" \
  "$P"
cat outputs/STAGE05E3F_GATE.json
column -t -s $'\t' outputs/STAGE05E3F_FAMILY_EXPANSION_SUMMARY.tsv
column -t -s $'\t' outputs/STAGE05E3F_NONWITNESS_ARCHITECTURE_EXPANSION.tsv
column -t -s $'\t' outputs/STAGE05E3F_RT_CONTEXT_SUMMARY.tsv
cat outputs/STAGE05E3F_RT_RESOURCE_STATUS.json
cat outputs/STAGE05E3F_MANUSCRIPT_NOTES.md
./collect_stage05e3f.sh
```

Scientific boundary: no E1/E2/E3D family membership, Stage05B score/rank/threshold, or RT-HMM result is changed. The frozen ten-accession fallback is used only to prove zero same-genome overlap when the sets are disjoint.
