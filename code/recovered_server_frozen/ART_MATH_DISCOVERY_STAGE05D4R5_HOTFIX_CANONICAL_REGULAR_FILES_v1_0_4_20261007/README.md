# Stage05D4-R5 v1.0.4 canonical regular-file layout hotfix

Purpose: fix a downstream compatibility bug left by v1.0.3. The v1.0.3 canonical query directories were real directories, but the files inside them were symlinks. Therefore `find ... -type f` returned 0 even though `find -L` returned 43. This hotfix converts each file-level symlink into a true regular file at the canonical path, preferring hard links (zero extra data copy) and falling back to a normal/reflink copy.

This hotfix does **not** rerun Boltz, alter P1/P2, change MSA/model/seed/sampling, recompute DNA scores/ranks, or filter candidates.

Run from any directory:

```bash
materialize_canonical_prediction_files.sh /path/to/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007
```

Expected gate: `PASS_BOLTZ2_P1_43_CANONICAL_REGULAR_FILE_LAYOUT` with 43 prediction dirs, 43 regular CIF files, and 43 regular confidence JSON files.
