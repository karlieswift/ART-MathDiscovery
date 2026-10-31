# Stage05D4-R5D v1.0.2 postprocessing hotfix

Scope: implementation-only/postprocessing-only repair after the successful P2 clean-RT Foldseek search produced `raw/P2_CLEAN_RT.m8` (65 rows) but the original summarizer stopped with `NameError: name 'Path' is not defined`.

Change: add `from pathlib import Path` to the summarizer. The hotfix also adds strict checks that the 40 confidence JSON query IDs exactly equal the frozen 40 P2 query IDs and that every Foldseek query maps back to a frozen P2 query.

No Boltz prediction rerun. No Foldseek rerun. No threshold/reference/candidate/DNA change.
