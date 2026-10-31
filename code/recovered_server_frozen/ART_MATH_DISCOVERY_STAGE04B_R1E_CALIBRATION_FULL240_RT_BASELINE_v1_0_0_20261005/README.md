# Stage04B-R1E: reference capability and full240 RT-only baseline

This is a new post-result analysis, NOT a replacement for R1D2 and NOT a new-genome discovery run. It preserves the 240 DNA genome universe and all 383 mathematical candidate regions. It neither modifies nor reads UG27 outputs.

## Purpose

Separate the team's DNA-region prioritization (ACNV + positional recurrence) from Pfam/HMMER homology evidence. Assess whether the current RT profiles detect the registered ART-associated RT references. Establish the missing whole-genome RT-only comparator over all240 genomes, including66 with no dual-arm candidates.

## Run on as-4125

Use the existing `art_rt_evidence` conda environment (Python3.11, HMMER3.4, Prodigal2.6.3). No new dependencies and no downloads. The two tool binary SHA256 values must also match the returned R1D2 run; the existing as-4125 environment is the intended environment.

```
conda activate art_rt_evidence
./run_selftest.sh
./run_stage04b_r1e.sh /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1 4
```

Source FASTAs are located ONLY using the exact frozen Stage04A manifest. Each file's byte SHA256 is checked. The earlier174 Prodigal results are copied into this separate directory after matching522 frozen raw-file hashes. Only the other66 genomes receive new Prodigal predictions. All240 predicted proteomes are scanned, irrespective of candidate or RT annotations.

Existing official neighborhood proteins are included as a separately labeled auxiliary branch, not treated as uniformly complete official annotations for all240 genomes. The primary full-universe comparator is predictor-only, using the same Prodigal implementation on every genome. It does NOT estimate true RT recall or ART sensitivity.

All targets (full240 prediction + existing official neighborhood proteins + six RT references) are exact-AA deduplicated in one frozen search universe. Three HMM profiles retain the exact R1D2 bytes, thresholds and seed. HMMER `-Z 100000` is fixed for this new comparison. E-values must NOT be reported as byte-exact replays of the smaller R1D2 search. Primary GA bit thresholds are unchanged.

## Reference calibration

The six reference proteins derive from previously registered coding coordinates and frozen genomes: Madawaska, Machias, SA1, Sen2, Ec78 and DRT2. Their nucleotide sequence, translation convention and AA hashes are bundled. No internal stops are replaced and no candidate-dependent translation is selected.

All six are RT references. The non-array RT systems are NOT RT-negative controls. This six-case panel measures limited detectability and system-type contrasts; it does not establish RT specificity or independent catalytic function. Madawaska/Machias/SA1 are related and were studied previously. MarsHill, L0050 and UG27 are not silently added to this reference denominator.

The unchanged R1D2 exact repeat algorithm is applied to three registered ART repeat-core contexts. A core is recovered only when a high-complexity group has at least three actual copies overlapping the registered core; a span crossing the core alone is insufficient. This is truth-informed evaluation of a frozen repeat assay, not truth-informed mathematical candidate scoring.

A reference not detected is a model/assay coverage finding, not a technical execution failure. No HMM thresholds, repeat thresholds or reference labels are adjusted after scores.

## Outputs

- `outputs/REFERENCE_RT_CAPABILITY.tsv`: six RT reference assay outcomes.
- `outputs/REFERENCE_REPEAT_CALIBRATION.tsv`: three repeat-core outcomes.
- `outputs/FULL240_PREDICTED_RT_ASSESSMENT.tsv`: all predicted genes and fixed profile evidence.
- `outputs/CORE_RT_LOCI_INSIDE_AND_OUTSIDE.tsv`: core-profile-positive loci inside/outside candidate neighborhoods, including zero-candidate genomes.
- `outputs/ALL383_BASELINE_LINKAGE.tsv`: all383 original candidates with new links; original fields/order retained.
- `outputs/COVERAGE_SUMMARY.json`: actual merged DNA footprints, not region-count compression.
- `outputs/POSITION_NULL_SUMMARY.json`: exploratory same-genome/same-width randomized-position comparison.
- `GATE_STAGE04B_R1E.json`: technical completion, scope limits, hashes. No automatic superiority/ART claim.

The position null uses1000 deterministic replicates and seed20261005. It preserves each region's length and genome assignment, but NOT GC, coding fraction, spacing, overlap configuration or merged footprint. It cannot be interpreted as a fully matched-workload comparison. It is conditional on the already selected genome distribution. Outcomes were previously inspected in this project; any p-value is exploratory, not untouched confirmation.

## Resume and collection

Run the SAME command to resume. Valid completed prediction/HMM outputs are reused. Changed source bytes, tool binaries or execution contract cause HOLD, not silent reuse. Do not modify the old results or delete previous HOLDs. This package uses its own lock.

After PASS or HOLD:

```
./collect_results.sh
```

Send the generated RESULTS archive and SHA256. It includes new66 prediction outputs; the already returned174 raw predictions are not duplicated.

`run_precompute_only.sh` runs only DNA-footprint accounting and known-core repeat calibration, requiring standard Python only. It does not run HMMER or Prodigal. The local validation file explicitly separates executed computations from pending server assays.
