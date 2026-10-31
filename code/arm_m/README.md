# ART Mathematical Genome Discovery — Stage03A frozen-engine L3 validation

This package freezes the final MarsHill-development engine and then performs **sequence-only, genome-blind validation** on the three Stage00-predeclared genomes not used in Stage02 development:

- MW349129 / Madawaska
- MW349128 / Machias
- MW218148.1 / SA1

No ART truth is read anywhere in this package. There is intentionally no truth-diagnostic script.

## Requirements
Python environment with numpy, pandas and the dependencies already used by the frozen GLMY/irregular engines. `bioh200` is the recommended environment. Internet access is needed only if the three FASTAs are not already placed under `inputs/`.

## Selftest
```bash
./run_selftest.sh
```

## Run
```bash
./run_stage03a_blind_validation.sh 32 0
```

Arguments: `32` workers, `0` removes large raw primary feature matrices after each genome is frozen.

If NCBI fetch is blocked, manually place raw single-record FASTAs as:
```text
inputs/MW349129.frozen.fasta
inputs/MW349128.frozen.fasta
inputs/MW218148.1.frozen.fasta
```
and rerun. They must be A/C/G/T-only and at least 4096 bp. The engine deliberately fails instead of inventing a new N-handling policy.

## Expected final gate
`PASS_STAGE03A_PREDECLARED_L3_BLIND_VALIDATION_SCORES_FROZEN`

Only after that gate is returned and audited may Stage03B reveal/import post-cutoff ART truth for evaluation.
