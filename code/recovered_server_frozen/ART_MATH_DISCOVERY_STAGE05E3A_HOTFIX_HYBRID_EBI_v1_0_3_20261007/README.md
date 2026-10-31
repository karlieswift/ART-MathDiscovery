# Stage05E3A hybrid EBI retrieval hotfix v1.0.3

Purpose: recover the 17 accessions that remained missing after v1.0.2, without changing the frozen E2A family/candidate universe or downstream biological analysis.

Observed v1.0.2 state: 44/61 records validated through ENA. The 17 missing IDs comprise sixteen `NC_*` RefSeq accessions plus `OK905446.1`.

Transport policy:
- Preserve all already validated v1.0.2 `.gbff` records.
- Route RefSeq nucleotide accessions to EMBL-EBI Dbfetch database `refseqn`, format `refseq`, raw style.
- Retry non-RefSeq missing accessions through EMBL-EBI ENA Dbfetch using both query and REST forms, versioned then versionless.
- Every recovered record must parse as a single non-empty Biopython GenBank record before entering `raw/genbank`.
- After 61/61 validation, call the original Stage05E3A `scripts/analyze_official_context.py` unchanged.

Scientific settings changed: **false**.

The original v1.0.2 failure is expected for `NC_*` accessions because `ena_sequence` is the INSDC primary nucleotide database, while EMBL-EBI exposes NCBI RefSeq nucleotide records separately via `refseqn`.
