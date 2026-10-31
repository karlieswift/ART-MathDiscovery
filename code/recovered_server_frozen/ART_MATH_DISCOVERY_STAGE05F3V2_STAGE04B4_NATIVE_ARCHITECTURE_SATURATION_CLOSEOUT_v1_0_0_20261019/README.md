# ART MathDiscovery Stage05F3V2 Stage04B4

Native-architecture + saturation closeout for the independent 54.38-Gb delta universe.

## What is frozen before looking at native architecture

- 54.377819157 Gb / 254 new accessions, with all original 131 source accessions excluded.
- Old/delta combined family partition and delta-only family space.
- Historical four-axis priority rule -> exact 32 delta families.
- Reference-set structural compatibility definition (2026-10-19): 3 STRONG + 9 SUGGESTIVE = 12 eligible families.
- Native context, repeat, strict ART-geometry and partner rules in `config/CONTRACT.json`.

## Run

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"

tar -xzf ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4_NATIVE_ARCHITECTURE_SATURATION_CLOSEOUT_v1_0_0_20261019.tar.gz
cd ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4_NATIVE_ARCHITECTURE_SATURATION_CLOSEOUT_v1_0_0_20261019

./run_selftest.sh
bash ./run_stage04b4_closeout.sh "$P" 32
```

The primary run is fail-closed and atomically commits `outputs/` only after a full PASS. It requires discoverable family membership for all 12 frozen structural-positive families. Every member candidate is replayed from its frozen raw-DNA contig using the frozen candidate coordinates and checked for exact translated length/no internal stop; the 12 frozen panel representatives are additionally required to match the pre-existing Stage04B2/B3 amino-acid FASTA byte-for-sequence. If any membership, source contig, coordinate replay, or representative sequence is unresolved, the run stops instead of silently falling back to representative-only analysis.

A representative-only technical diagnostic is available only by directly invoking the Python driver with `--representative-only`; it is not the primary closeout.

## Expected primary output

- `outputs/DELTA12_STRUCTURAL_POSITIVES.tsv`
- `outputs/DELTA12_FAMILY_MEMBERSHIP.tsv`
- `outputs/DELTA12_NATIVE_CONTEXTS.tsv`
- `outputs/DELTA12_NATIVE_CONTEXTS.fna`
- `outputs/DELTA12_REPEAT_ARRAYS.tsv`
- `outputs/DELTA12_NEAREST_DOWNSTREAM_ORFS_WITH_CLUSTER.tsv`
- `outputs/DELTA12_FAMILY_NATIVE_ARCHITECTURE_SUMMARY.tsv`
- `outputs/DELTA12_STRATUM_ARCHITECTURE_SUMMARY.tsv`
- `outputs/STAGE04B4_CLOSEOUT_GATE.json`
- `outputs/STAGE04B4_OUTPUT_SHA256.tsv`

Then:

```bash
./collect_results.sh
```

Upload the resulting `ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4_RESULTS_*.tar.gz` and sidecar SHA256 for audit/manuscript integration.
