# Stage05D4-R5 v1.0.3 canonical prediction-layout hotfix

Purpose: repair a downstream path/layout incompatibility only.

The v1.0.2 collector correctly assembled 43/43 Boltz predictions but represented each
`boltz_out/predictions/<query>` as a symlink to a directory. The frozen downstream
`run_full3d_foldseek_and_summarize.sh` used GNU `find` without `-L`, so it did not
traverse those command-line symlink directories and reported `expected 43 flat CIFs, got 0`.

This hotfix creates real canonical query directories under `boltz_out/predictions/` and
places file-level symlinks to the already completed Boltz outputs inside them. It preserves
the prior v1.0.2 assembly and does not rerun Boltz.

Scientific settings/candidate universe/DNA scores/ranks are unchanged.
