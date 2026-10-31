# Stage05D4-R5 v1.0.2 collector-only hotfix

The four independent single-GPU Boltz workers completed successfully, but v1.0.1 assumed outputs lived directly at `boltz_out_shards/gpuN/predictions`.
Boltz 2.2.x may nest results under an additional `boltz_results_*` directory.

This hotfix recursively locates exactly one `predictions/` directory beneath each shard, assembles canonical `boltz_out/predictions`, and validates 43 CIF + 43 confidence JSON files.
It does not rerun Boltz, change P1/P2, change any model or sampling setting, or modify DNA scores/ranks.
