# ART Mathematical Genome Discovery — Stage04B-R1 post-freeze biological annotation

This stage starts the biological unblinding of the **383 frozen Stage04A-R1 ACNV + generic positional recurrence candidate regions**.

It does **not** recompute any mathematical score, rank, q-value, lag, window, merge, or threshold. It does **not** discard candidates. The complete 383-row registry is carried forward.

## What this stage does

1. Verifies the frozen Stage04A universe manifest and the exact Stage04A-R1 383-candidate registry.
2. Reads the exact frozen Stage04A FASTA bytes from the original Stage04A project root.
3. Downloads NCBI GenBank records for the 174 versioned accessions represented by the 383 candidates.
4. Requires the GenBank sequence to match the Stage04A frozen FASTA SHA256 before using its annotation.
5. Extracts every candidate sequence and a fixed +/-10 kb biological neighborhood.
6. Exports CDS/RNA/repeat annotations around every candidate.
7. Performs a **descriptive RT-keyword triage only**. No candidate is removed if it lacks an RT keyword.
8. Exports all neighborhood proteins for the next sequence-based RT/domain scan, which is required because novel RTs can be annotated as hypothetical proteins.

## Why this is the correct next step

Stage04A-R1 found 383 frozen two-arm mathematical candidates in 174/240 sequences. They are candidate regions, not ART identity calls. This stage asks whether those frozen candidates have biological neighborhoods compatible with reverse-transcriptase-associated array systems while preserving the pre-annotation mathematical freeze.

## Run

```bash
conda activate bioh200
cd ART_MATH_DISCOVERY_STAGE04B_R1_POSTFREEZE_BIOLOGICAL_ANNOTATION_v1_0_0_20261004
./run_selftest.sh

./run_stage04b_r1.sh \
  /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE04A_L4_BROAD_UNIVERSE_FROZEN_DISCOVERY_v1_0_0_20260928
```

Optional NCBI API key:

```bash
export NCBI_API_KEY='...'
```

## Main outputs

- `outputs/STAGE04B_R1_SEQUENCE_PROVENANCE.tsv`
- `outputs/STAGE04B_R1_CANDIDATE_NEIGHBORHOODS.tsv`
- `outputs/STAGE04B_R1_NEIGHBORHOOD_FEATURES.tsv`
- `outputs/STAGE04B_R1_CANDIDATE_BIOLOGY_SUMMARY.tsv`
- `outputs/STAGE04B_R1_RT_KEYWORD_TRIAGE.tsv`
- `outputs/STAGE04B_R1_CANDIDATE_REGIONS.fasta`
- `outputs/STAGE04B_R1_NEIGHBORHOODS.fasta`
- `outputs/STAGE04B_R1_NEIGHBORHOOD_PROTEINS.faa`
- `outputs/STAGE04B_R1_SUMMARY.json`
- `outputs/STAGE04B_R1_GATE.json`

Expected clean gate:

`PASS_STAGE04B_R1_ALL_383_FROZEN_CANDIDATES_ANNOTATED_WITH_EXACT_SEQUENCE_PROVENANCE`

If any NCBI GenBank sequence does not match the Stage04A frozen FASTA bytes, the stage fails closed and biological interpretation is not authorized for that run.
