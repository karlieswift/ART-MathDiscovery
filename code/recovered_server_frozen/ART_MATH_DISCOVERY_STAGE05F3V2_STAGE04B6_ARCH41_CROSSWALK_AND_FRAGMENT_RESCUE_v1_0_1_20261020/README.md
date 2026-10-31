# Stage05F3V2 Stage04B6 — Architecture41 crosswalk and fragment rescue

This stage acts on the frozen 41 Stage04B5 strict+recurrent-partner families.

- 23 MIXED families: direct historical 6.12-Gb family-space crosswalk with the same 30% identity / 70% coverage rule.
- 18 DELTA_ONLY families: deterministic longest-member audit, motif-C-like descriptive audit, Logan-unitig context recovery, and same-frame maximal stop-to-stop fragment extension.
- No Boltz/Foldseek is run here. Structure is deferred until fragment recovery is closed.

Fast path:

```bash
./run_stage04b6_prepare.sh "$P" 32
./run_stage04b6_unitig_rescue.sh "$P" 4
./collect_results.sh
```

The unitig runner downloads only the unique accessions selected by the frozen longest-member rule (expected 11 unique accessions for the 18 delta-only families).
