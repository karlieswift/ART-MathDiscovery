# ART MathDiscovery Stage05E3A — E2A official context + lineage audit

Purpose: perform post-discovery biological-context audit of the 12 frozen Stage05E2 E2A exact-primary joint-signature architecture families.

Frozen before GenBank retrieval:
- 12 families
- 119 candidate loci
- 61 sequence accessions
- 80 candidate-level family-exclusive exact-primary repeat + compact-FAM50 witness pairs spanning 44 genomes
- No structure result is used for selection; no DNA score/rank is recomputed; no candidate is filtered.

The package fetches official NCBI GenBank records, records taxonomy/source metadata, audits ±5 kb and ±10 kb official CDS context, maps E2 witness Prodigal ORFs to official CDS by coordinate/strand overlap, and flags explicit curated retro-system annotation keywords. Absence of a keyword is not novelty evidence. Presence is only context/triage evidence.

Run:
```bash
conda activate art_novelty
./run_selftest.sh
./run_stage05e3a.sh
cat outputs/STAGE05E3A_SUMMARY.json
cat outputs/GATE_STAGE05E3A.json
column -t -s $'\t' outputs/STAGE05E3A_FAMILY_BIOLOGICAL_CONTEXT.tsv | head -30
./collect_stage05e3a.sh
```

Optional NCBI environment variables: `NCBI_API_KEY`, `NCBI_EMAIL`.
