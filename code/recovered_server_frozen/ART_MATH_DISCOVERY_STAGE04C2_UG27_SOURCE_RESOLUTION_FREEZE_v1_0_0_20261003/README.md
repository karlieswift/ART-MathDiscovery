# ART MathDiscovery Stage04C2 — UG27 exact experimental source-resolution freeze

Purpose: resolve the **exact five experimentally tested UG27 systems** from the 2026 Minerva study **before any ACNV/recurrence scoring**.

This package is provenance-only. It MUST NOT run the frozen detector, retune lag/q thresholds, or substitute arbitrary UG27 representatives.

## Expected workflow

1. Obtain the official Minerva supplementary table from Zenodo record `22903368`.
2. Place it at `inputs/Minerva_SupplementaryTables_preprint.xlsx` (or run the download helper).
3. Run self-test.
4. Run source-resolution scan.
5. Inspect `outputs/GATE.json` and the candidate tables.

The gate is deliberately conservative:
- `PASS_STAGE04C2_UG27_EXPERIMENTAL_FIVE_EXACT_SOURCE_IDS_FROZEN` only when exactly five experimental systems are uniquely resolved to exact source identifiers.
- Otherwise `HOLD_STAGE04C2_UG27_EXACT_SOURCE_IDS_UNRESOLVED`.

A HOLD is scientifically correct and does **not** authorize scoring.
