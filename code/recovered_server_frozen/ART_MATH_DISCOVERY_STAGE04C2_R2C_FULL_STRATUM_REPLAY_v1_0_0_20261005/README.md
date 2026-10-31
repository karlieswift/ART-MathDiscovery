# Stage04C2-R2C: complete-stratum, resource-bounded replay

Version 1.0.0 / 2026-10-05. Linux, Python >=3.10, NumPy and pandas; use bioh200.
This package does not install or upgrade dependencies and does not edit Stage04A.

## What changes, and what does not

The original mathematical functions are imported from an isolated snapshot of
the SHA256-locked Stage03A ZIP. Stage04A's runtime source files must also match
that ZIP. Extraction is scheduled by COMPLETE (window length, strand) stratum.
Each full stratum is scored with the original robust_score; ranks, strand
collapse, empirical calibration, recurrence and region merge use the complete
original populations. A 256-item processing batch is only a task-submission
batch, NEVER a calibration population. Exact median is not approximated.

All original computational branches are retained for reproducibility; the
reported two-arm view is ACNV + generic positional recurrence. Original
historical region ranks are retained and are not re-labelled as two-arm ranks.
This implementation is not the original entrypoint byte-for-byte: it is an
implementation variant that MUST pass the comparison below before UG27.

All 5 assemblies / 15 contigs / 20,664,383 bp / 2,542,244 oriented windows stay
in scope. Identical contigs in different assemblies are NOT silently removed.
No contigs are concatenated. No locus coordinates, plasmid alignments or truth
annotations are read by the scoring runner.

## Run

Set the two exact existing source roots; no glob-based root selection is used.

```bash
conda activate bioh200
A04=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE04A_L4_BROAD_UNIVERSE_FROZEN_DISCOVERY_v1_0_0_20260928
UG=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE04C2_UG27_SOURCE_RESOLUTION_FREEZE_v1_0_1_20261003
./run_selftest.sh
./run_resource_plan.sh "$A04" "$UG"
./run_historical_validation.sh "$A04" 8
./run_ug27_score.sh "$A04" "$UG" 8
```

The score command refuses to start without a current PASS validation gate. It
checks source, implementation, runtime and validated-board hashes again.
Historical validation accessions are frozen in config: NC_002666.1 and
NC_028671.2. They are not chosen by looking at UG27 scores. Three full TSV boards
per accession must have EXACT serialized content equality (gzip decompressed
for that comparison). Compressed container equality is also reported separately.
No numeric tolerance can silently turn a difference into a PASS.

Validation and scoring write detailed progress to logs/validate_ACCESSION.log
and logs/ug27_SCORING_ID.log. Re-running the same score command resumes completed
complete-stratum score checkpoints and completed contigs after verifying hashes.
A failed/incomplete group is recomputed; completed inputs and historical boards
are never deleted. If an interrupted process leaves .run_lock, verify no process
is still active before removing that empty lock directory.

## Resource interpretation

91,907,205,088 bytes is the cumulative size of the two raw arrays across ALL 15
contigs, not a compulsory simultaneous RAM or disk requirement. Sequential
original execution would peak at 24,953,411,872 raw-array bytes on the largest
contig. Here the largest single extraction memmap is 3,834,227,216 bytes.
These are ARRAY STORAGE estimates, not total resident RAM: robust_score still
allocates intermediate matrices, and ranking/recurrence have additional costs.
All feature extractions remain required, so this package does not promise a
proportional CPU speedup or a completion time. Contigs run sequentially; inner
process count defaults to 8. No GPU implementation or precision change is made.

## Outputs and return package

- outputs/HISTORICAL_EQUIVALENCE_GATE.json
- outputs/RESOURCE_PLAN.json and RESOURCE_PLAN.tsv
- outputs/contigs/SCORING_ID/ (all original-format region and window boards)
- outputs/TWO_ARM_REGION_VIEW.tsv (all dual-supported regions, no Top-K)
- outputs/UG27_SCORE_FREEZE_GATE.json
- outputs/SCORE_OUTPUT_MANIFEST.tsv

```bash
./collect_results.sh
```

Return the generated RESULTS tar.gz and its .sha256. Large temporary feature
matrices, the source archive and the Stage04A original runtime are excluded.
Collect only after a run exits; the collector refuses an active run lock.

## Test status and scientific boundary

Local package tests use synthetic data and the displayed mathematical functions,
including a multi-process extraction/replay and checkpoint-resume fixture. They
are NOT the real historical-conformance test and NOT a completed UG27 run.
Both real historical genomes must pass on the user's machine. Full UG27 scores
must then be collected and hashed before a separate post-score evaluation.

This is source-known, score-blind frozen-parameter replay, not a claim that the
analyst has never seen the loci. R3 supplied cassette correspondence: four exact
long alignments and one near-exact alignment. A plasmid-insert/cassette alignment
is not, by itself, an independently defined ncRNA-array truth interval.
