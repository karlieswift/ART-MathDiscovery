# Stage05E3A v1.0.2 — ENA provider hotfix

Purpose: recover the frozen Stage05E3A official-context analysis after the server IP was redirected by NCBI to `misuse.ncbi.nlm.nih.gov/error/abuse.shtml` for both EFetch and Sequence Viewer.

This hotfix changes only the retrieval provider/transport layer. It does **not** change the frozen E2A families, 119 candidate loci, 61 sequence accessions, witness definitions, candidate membership, DNA scores/ranks, biological keyword rules, or downstream analysis logic.

The script retrieves the same INSDC nucleotide records from EMBL-EBI ENA using EBI Dbfetch (`ena_sequence`, raw EMBL flatfile), validates them with Biopython, converts them deterministically to GenBank flatfiles, and then invokes the unchanged Stage05E3A `analyze_official_context.py`.

Primary route:
`https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ena_sequence&id=<accession>&format=embl&style=raw`

Fallback: accession without version through the same Dbfetch endpoint.

The original NCBI-blocked payloads, if present, are preserved under `raw/genbank_ncbi_blocked_v100_v101/` before new retrieval begins.

Run:

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
E3A="$P/ART_MATH_DISCOVERY_STAGE05E3A_E2A_OFFICIAL_CONTEXT_LINEAGE_v1_0_0_20261007"
"$P/ART_MATH_DISCOVERY_STAGE05E3A_HOTFIX_ENA_PROVIDER_v1_0_2_20261007/run_stage05e3a_ena_v102.sh" "$E3A"
```
