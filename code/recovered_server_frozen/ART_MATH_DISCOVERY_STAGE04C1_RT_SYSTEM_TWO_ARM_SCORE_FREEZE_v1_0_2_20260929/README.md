# Stage04C1 RT-system two-arm score freeze v1.0.2

Hotfix over v1.0.1. The C0 registry column is `whole_sequence_accession`; v1.0.1 incorrectly referenced the nonexistent `sequence_accession` field and therefore failed during preflight before any scoring. v1.0.2 fixes only that preflight schema/provenance lookup. Frozen ACNV source, recurrence definition, thresholds, eligible controls, R3 sequence-freeze requirement, and scoring code are unchanged.

Scientific status: no control scoring occurred in the failed v1.0.1 run.
