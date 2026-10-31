# Stage04B7 — rescued-delta structural closure + mixed23 family-wide replication closure

Expected starting point: Stage04B6 v1.0.1 is complete with 12 `UNIQUE_STOP_BOUNDED_BOTH`, 1 `UNIQUE_EDGE_CENSORED`, and 5 `ABSENT` delta-only families.

Run order:

1. `./run_stage04b7_prepare.sh "$P" 32`
2. `conda activate art_structure && ./run_stage04b7_boltz_shards.sh "$P" 0,1,2,3`
3. `conda activate art_novelty && ./run_stage04b7_foldseek_finalize.sh "$P" 32`
4. `./collect_results.sh`

The Boltz step uses independent single-GPU shards, not multi-GPU DDP.
