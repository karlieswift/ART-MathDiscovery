# Stage04B7R1 technical resume + replication semantic guard

Purpose:
1. Resume only missing Stage04B7 Boltz2 queries with the exact frozen structural protocol; do not rerun successful queries.
2. Audit the MIXED23 historical5 crosswalk for many-to-many fanout. The existing 30% identity / 70% coverage pairwise search is valid as sequence-space compatibility evidence, but many-to-many mapping does not by itself establish one-to-one historical system replication.

No science thresholds are changed. No candidate membership is changed. No retuning is performed.

Usage:

```bash
bash scripts/run_b7r1_resume_missing_boltz.sh <B7_ROOT> <PROJECT_ROOT> 0,1,2,3
python scripts/audit_b7_replication_semantics.py <B7_ROOT>
```
