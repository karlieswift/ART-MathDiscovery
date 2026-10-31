# Stage05D4-R5D tmp-directory hotfix v1.0.1

Observed failure: Foldseek `easy-search` stopped with `Cannot create temporary directory tmp/p2_clean_rt` after the script had already validated the informative R5C calibration and staged all 40 P2 CIF files.

Root cause: `run_clean_rt_screen.sh` removes `tmp/p2_clean_rt` but creates only `work/p2_cif_flat_regular`, `raw`, and `logs`; it does not ensure that the parent `tmp/` directory exists. Foldseek cannot create the nested path if its parent is absent.

This wrapper creates/preserves a usable `tmp/` parent, removes only the R5D Foldseek scratch subdirectory, and then invokes the original frozen R5D screen script unchanged.

Scientific settings changed: no. Boltz predictions rerun: no. P2 universe changed: no. R5C clean-RT reference changed: no. Foldseek thresholds changed: no. DNA scores/ranks changed: no.
