# ART Mathematical Genome Discovery — Stage04A-R1

## Purpose

Stage04A-R1 is a **post-freeze deterministic two-arm projection** of the completed Stage04A broad-universe region board. It does not rescore sequences, recompute q values, recompute region ranks, change thresholds, use biological annotation, use truth labels, or authorize retuning.

The paper-facing two-arm definition is fixed as **ACNV + generic positional recurrence**. A region is retained when both `ACNV_FAMILY` and `GENERIC_RECURRENCE` are already present in the frozen region board's `support_families` field.

## Canonical source locks

- `STAGE04A_ALL_FROZEN_UNION_REGIONS.tsv` SHA256: `272b59fd09d853bc7336e364a64e909c3cca7b082a684c69355d873707be9337`
- `STAGE04A_UNIVERSE_MANIFEST.tsv` SHA256: `24f69c4f463678413d772cc567332bd0b94ff18019d53ea7e454aa5e30fb5ea4`

## Frozen closeout result

- 240 sequences
- 5,838 frozen union regions
- 383 ACNV+recurrence regions
- 6.5605% of the frozen region universe
- 15.2428-fold candidate-space compression
- 174/240 sequences contain at least one two-arm region
- 66/240 contain none
- pooled two-arm burden: 23.1074 regions/Mb
- median frozen within-sequence rank among two-arm regions: 2 (IQR 1–4)
- 142 sequences have a two-arm region at frozen rank 1

These 383 regions are **mathematical candidates for downstream biological annotation**, not ART identity calls. No global ART probability or post-hoc rank threshold is defined here.

## Outputs

- `outputs/STAGE04A_R1_DUAL_ARM_REGIONS.tsv`: all 383 retained regions; this is the primary Stage04A-R1 candidate registry.
- `outputs/STAGE04A_R1_PER_SEQUENCE_DUAL_ARM_BURDEN.tsv`: per-sequence burden summary.
- `outputs/STAGE04A_R1_BEST_DUAL_ARM_PER_SEQUENCE.tsv`: one descriptive best frozen-rank two-arm region per positive sequence. This is not a survival cutoff; all 383 regions remain primary.
- `outputs/STAGE04A_R1_SUMMARY.json`
- `outputs/STAGE04A_R1_GATE.json`

## Run

```bash
./run_selftest.sh
./run_stage04a_r1.sh /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE04A_L4_BROAD_UNIVERSE_FROZEN_DISCOVERY_v1_0_0_20260928
```

Expected formal gate:

`PASS_STAGE04A_R1_FROZEN_TWO_ARM_UNIVERSE_CLOSEOUT`

## Scientific boundary

Stage04A-R1 establishes a frozen mathematical candidate set. Determining whether any candidate is ART-like, RT-associated, array-associated, previously known, or novel belongs to the next biological-annotation stage. Stage04A-R1 itself makes no such identity claim.
