# ART-MathDiscovery

Code-only repository for the ART-MathDiscovery computational workflows.

## Contents

- `code/arm_m/`: ACNV/recurrence engine and validation code.
- `code/clean_two_arm_replay/`: clean two-arm replay code.
- `code/formal_baseline/`: baseline inventory and scoring code.
- `code/ug27_preunblind/` and `code/ug27_exact_truth/`: UG27 analysis code.
- `code/j1_j2_supplement/`: J1/J2 calibration, cohort, and test code.
- `code/recovered_server_frozen/`: frozen historical source-code snapshots.
- `config/`: authoritative runtime configuration.
- `scripts/`: repository checks and source self-tests.

Paper sources, manuscript files, execution logs, generated outputs, provenance
records, and internal process materials are intentionally excluded.

## Quick check

```bash
python3 scripts/verify_release.py
python3 -m pip install -r requirements-tests.txt
bash scripts/run_selftests.sh
```

Some full workflows require external datasets that are not distributed in this
code-only repository.
