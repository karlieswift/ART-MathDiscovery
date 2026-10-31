Stage05F3B4E G1-R2 HTTPS Common URL API transport hotfix

Reason:
- BLAST+ -remote failed on as-4125 and r750-1 with the same CRPCClient Blast4-request bad-stream error.
- Failure happened before any G1 search result, so transport may be replaced result-blind.

Scientific search scope is unchanged:
- RT5: blastp nr, E<=1e-3, 200 targets, low-complexity filtering on.
- partner5: same.
- architecture5: megablast nt, E<=1e-10, 200 targets, filtering on; NO_DATABASE_OVERRIDE=true preserves nt scope.
- Existing summarize_global.py and novelty thresholds are byte-identical to base F3B4E.

This hotfix uses the official NCBI Common URL API over HTTPS and retrieves XML2, then converts XML2 to the original frozen TSV columns.
