# Stage04C2-R4 exact UG27 ncRNA-array truth resolution

Purpose: final post-hoc truth refinement for the five prespecified hCom2 UG27 systems. This stage does **not** retune, recompute, or reinterpret the frozen detector board. It first attempts to freeze exact author-defined ncRNA units independently of ACNV/recurrence results; only if that truth gate passes may a separate R4B evaluation read the frozen R2D clean board.

## Why this is a truth-resolution stage, not a detector stage

The Minerva paper defines UG27 ncRNA units using two iteratively refined covariance models. Published Methods state an E-value threshold of 0.01 for ncRNA units and resolve overlaps between the two models by retaining the lower-E-value match. The public Zenodo publication record exposes the joined manuscript PDF and SupplementaryTables workbook; the public Minerva GitHub repository is framework code. As of package freeze, a machine-readable pair of the final UG27 covariance models or a five-system native coordinate table had not been independently verified in those public materials. Therefore this package refuses to manufacture 'exact' boundaries from motifs, generic intergenic regions, or our detector output.

## Accepted exact-truth sources

1. `--truth-tsv`: author machine-readable exact unit coordinates with columns: `system source_contig_id start0 end0 strand unit_id` (optional `evalue`, `source_method`).
2. `--cm-main` and `--cm-alt`: the two **author final UG27 covariance models**. The resolver runs Infernal `cmsearch`, uses published unit E<=0.01, resolves overlapping model matches by lower E-value, and keeps author-model units inside the frozen native cassette.
3. `--unit-fasta`: exact author-defined ncRNA unit sequences. Each unit must map uniquely within the corresponding frozen native cassette.

If none of these is available, R4A returns `HOLD_STAGE04C2_R4A_EXACT_ARRAY_TRUTH_NOT_RESOLVED`. That HOLD is scientifically preferable to a post-hoc guessed array boundary.

## Primary exact-array endpoint frozen before exact truth is inspected

R4B primary endpoint: strict 0-based half-open overlap between a frozen clean ACNV+GENERIC_RECURRENCE dual-arm region and **at least one exact author-defined ncRNA unit**. Array-span overlap, number/fraction of units overlapped, ACNV-only unit overlap, and nearest dual-arm distance are secondary/descriptive. No nearby threshold is defined.

## Run

```bash
conda activate bioh200
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
UG="$P/ART_MATH_DISCOVERY_STAGE04C2_UG27_SOURCE_RESOLUTION_FREEZE_v1_0_1_20261003"
R3="$P/ART_MATH_DISCOVERY_STAGE04C2_R3_POSTFREEZE_CASSETTE_UNBLINDING_v1_0_0_20261018"
R2D="$P/ART_MATH_DISCOVERY_STAGE04C2_R2D_PREUNBLIND_CLEAN_TWO_ARM_PROJECTION_v1_0_0_20261018"

./run_selftest.sh

# First inventory only: expected to HOLD if no exact author truth artifact is present.
./run_stage04c2_r4a_truth_resolution.sh "$UG" "$R3"

# If you have the two final author CMs:
./run_stage04c2_r4a_truth_resolution.sh "$UG" "$R3" \
  --cm-main /path/to/author_UG27_main.cm \
  --cm-alt /path/to/author_UG27_alt.cm

# OR an exact author coordinate table:
./run_stage04c2_r4a_truth_resolution.sh "$UG" "$R3" \
  --truth-tsv /path/to/author_exact_units.tsv

# Only after R4A PASS:
./run_stage04c2_r4b_exact_array_eval.sh "$R2D"
./collect_results.sh
```

## Frozen boundaries

- no ACNV/recurrence score used to resolve truth
- no detector region used to resolve truth
- no motif-only heuristic accepted as exact truth
- no intergenic boundary accepted as exact truth
- no score/q/window/lag/merge/rank recomputation
- no post-hoc nearby rescue threshold
- `science_change=false`
