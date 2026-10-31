# Stage05E3A final-record DDBJ recovery hotfix v1.0.4

Observed state after v1.0.3: 60/61 frozen sequence records are already present and valid. The sole missing accession is `OK905446.1`.

This package is a transport-only recovery for that one accession:
- Preserve all 60 already validated `.gbff` records.
- Retrieve only `OK905446.1` from DDBJ getentry, which serves the DDBJ/EMBL/GenBank INSD flatfile database and supports accession.version queries.
- Require the response to start with `LOCUS`, parse as exactly one non-empty Biopython GenBank record, and match accession base `OK905446` before installation into `raw/genbank/OK905446.1.gbff`.
- Revalidate all frozen 61 records.
- Call the original Stage05E3A `scripts/analyze_official_context.py` unchanged.

No scientific universe, family membership, DNA score/rank, threshold, structure-based selection rule, or downstream interpretation rule is changed.
