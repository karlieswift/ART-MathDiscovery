# Stage05E3A GenBank retrieval hotfix v1.0.1

Scope: implementation/transport only. The frozen 12 E2A families, 119 candidate loci, 61 sequence accessions, 80 witness loci, DNA scores/ranks, and all biological decision rules are unchanged.

Observed v1.0.0 failure: all 61 HTTP retrievals returned payloads that failed GenBank-flatfile validation, so `analyze_official_context.py` correctly never ran and no Stage05E3A outputs/gate were produced.

v1.0.1 changes only the retrieval layer:
1. use NCBI EFetch `rettype=gb&retmode=text` (documented GenBank flat-file route),
2. retry accession without version,
3. fallback to NCBI Sequence Viewer GenBank text endpoint,
4. atomically accept only Biopython-parseable GenBank payloads,
5. preserve old invalid payloads under `raw/genbank_v100_invalid/`,
6. write route-level diagnostics including response first line/content-type/bytes/hash.

No candidate filtering or scientific retuning is performed.
