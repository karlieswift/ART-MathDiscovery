# ART Mathematical Genome Discovery — Stage04A accelerated L4 broad-universe frozen discovery

Purpose: move immediately from the completed Stage03B 3/3 external validation to a larger, unknown sequence universe **without any retuning**.

## What is frozen

The exact Stage03A engine ZIP is embedded at `provenance/ART_MATH_DISCOVERY_STAGE03A_FROZEN_ENGINE_L3_VALIDATION_v1_0_0_20260927.zip` and must have SHA256:

`902bddb7813a6aa04736c4a5dc4cd2439821c8a21efb1c931a98ad0a9a3bfe3c`

Stage04A never uses Stage03B truth files, the revealed 15-nt consensus, known repeat positions, or RT/partner positions for scoring.

## Fastest route now: curated pilot

If NCBI `datasets` CLI is installed:

```bash
cd ART_MATH_DISCOVERY_STAGE04A_L4_BROAD_UNIVERSE_FROZEN_DISCOVERY_v1_0_0_20260928
./run_selftest.sh

# deterministic 250-genome complete RefSeq Caudoviricetes pilot
STAGE04_SAMPLE_N=250 ./download_and_prepare_curated_ncbi.sh

./run_stage04a_primary.sh \
  work/ncbi_curated/prepared/STAGE04A_L4_CURATED_COHORT.tsv \
  4 8 0
```

`4 8 0` means four genomes concurrently, eight Stage03A worker processes per genome, and no large raw matrices retained. If RAM pressure is high, use `2 16 0`. GPU is not the limiting resource for these frozen Python feature extractors.

## Genuine prospective arm

The temporal boundary is fixed at **2026-09-29 00:00:00 UTC**. Do not call pre-boundary sequences prospective.

```bash
./download_and_prepare_prospective_ncbi.sh
./run_stage04a_primary.sh \
  work/ncbi_prospective_20260929/prepared/STAGE04A_L4_PROSPECTIVE_COHORT.tsv \
  4 8 0
```

If no records are available yet, the downloader exits cleanly; retrying later does not change the frozen policy.

## Custom / metavirome cohort

Prepare a TSV with columns:

`sequence_id  arm  source  release_date  fasta_path  cluster_id  notes`

Each FASTA must be a single A/C/G/T-only record >=4096 bp. For `L4_PROSPECTIVE`, `release_date` must be ISO `YYYY-MM-DD` and >= 2026-09-29. The package fails closed on N because Stage03A had no N-handling policy.

Run:

```bash
./run_stage04a_primary.sh /absolute/path/cohort.tsv 4 8 0
```

## Main outputs

- `outputs/STAGE04A_UNIVERSE_MANIFEST.tsv`
- `outputs/genomes/<id>/STAGE03A_UNION_REGIONS.tsv` — complete frozen board for every sequence
- `outputs/STAGE04A_ALL_FROZEN_UNION_REGIONS.tsv` — all regions across the universe
- `outputs/STAGE04A_PER_SEQUENCE_BURDEN.tsv`
- `outputs/STAGE04A_COHORT_SUMMARY.json`
- `outputs/STAGE04A_FINAL_GATE.json`

Expected final gate:

`PASS_STAGE04A_FROZEN_ENGINE_BROAD_UNIVERSE_PRIMARY_SCORES_FROZEN`

Only after archiving/hashing this output should Stage04B inspect genes, RTs, partner proteins, repeats, literature, or RNA data.
