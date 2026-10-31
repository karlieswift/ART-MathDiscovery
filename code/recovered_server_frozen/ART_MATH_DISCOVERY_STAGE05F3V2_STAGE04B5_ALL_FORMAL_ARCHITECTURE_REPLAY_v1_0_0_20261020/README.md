# Stage05F3V2-Stage04B5 — 54G all-formal-locus architecture-first replay

Purpose: replay the final 6.12-Gb discovery ordering on **all 88,487 frozen formal RT-like instances** from the 54.3778-Gb delta universe, rather than only the 32 protein-priority families or 12 structural-positive families.

## What is frozen

- 54.3778-Gb / 254-accession source panel: unchanged.
- formal RT membership: 88,487, unchanged; no HMM rerun.
- six-frame candidate coordinates and amino-acid sequences: replayed exactly.
- repeat detector: k=8,10,12,14,16; >=3 copies; start-to-start gap 60–360 bp; exact shared extension <=80 bp.
- strict ART geometry: upstream; 3–21 copies; every inferred spacer 120–220 bp.
- partner scan: nearest candidate-oriented downstream complete ORF, code 11 start/stop semantics, >=30 aa, both strands, <=10 kb.
- MMseqs family/partner clustering: 30% identity, 70% coverage, cov-mode 0, cluster-mode 0, sensitivity 7.5.

## Important design change relative to earlier Stage04B4

This is not a threshold change. It changes the *analysis order* to match the final 6.12-Gb discovery route:

```
88,487 formal loci
  -> all-locus context/evaluability
  -> all-locus frozen repeat + partner scan
  -> family collapse
  -> MIXED replication arm + DELTA_ONLY novelty arm
  -> architecture-selected families
  -> later structure
```

## Runtime

The expensive part is streaming the 254 frozen `.contigs.fa.zst` sources once. Default `workers=4` scans up to four accessions in parallel. Re-running is resumable: a completed accession is skipped after its `DONE.json` is validated.

## Run

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"

tar -xzf ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B5_ALL_FORMAL_ARCHITECTURE_REPLAY_v1_0_0_20261020.tar.gz
cd ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B5_ALL_FORMAL_ARCHITECTURE_REPLAY_v1_0_0_20261020

./run_selftest.sh

bash ./run_stage04b5.sh "$P" 32 4 \
  2>&1 | tee logs/run_stage04b5_outer_console.log
```

Arguments: `project_root`, `mmseqs_threads`, `accession_workers`.

## Expected final gates

```
PASS_STAGE04B5_PREFLIGHT_88487_FORMAL_FROZEN
PASS_STAGE04B5_ALL_FORMAL_ACCESSION_SCAN_COMPLETE
PASS_STAGE04B5_FAMILY_AND_PARTNER_CLUSTERING_COMPLETE
PASS_STAGE04B5_54G_MATCHED_ARCHITECTURE_REPLAY_COMPLETE
PASS_STAGE04B5_ATOMIC_OUTPUT_COMMIT
```

Then:

```bash
./collect_results.sh
```
