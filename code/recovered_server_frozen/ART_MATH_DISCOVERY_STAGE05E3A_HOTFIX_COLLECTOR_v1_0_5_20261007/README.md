# Stage05E3A v1.0.5 collector-only hotfix

Purpose: repair the Stage05E3A result collector after the retrieval transport chain moved from the original single NCBI manifest to the validated multi-provider ENA + EBI RefSeq + DDBJ provenance.

This hotfix does **not** rerun retrieval or biological analysis and does not change the frozen scientific universe, family membership, thresholds, DNA scores/ranks, or interpretation rules.

The original `collect_stage05e3a.sh` hard-coded the obsolete file `raw/STAGE05E3A_GENBANK_RETRIEVAL_MANIFEST.tsv`. After v1.0.4, the canonical retrieval ledger is `raw/STAGE05E3A_COMBINED_PROVIDER_MANIFEST_v104.tsv`; the older provider-specific manifests/diagnostics are retained as optional audit evidence.

Use:

```bash
./collect_stage05e3a_v105.sh /path/to/ART_MATH_DISCOVERY_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_v1_0_0_20261007
```
