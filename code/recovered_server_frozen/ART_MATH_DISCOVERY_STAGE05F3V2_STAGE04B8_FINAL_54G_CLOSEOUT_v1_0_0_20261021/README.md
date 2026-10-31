# Stage04B8 Final 54G Closeout v1.0.0

Lightweight CPU-only closeout. It reads completed Stage04B5/B6/B7 outputs and produces the final 54G waterfall, old-vs-delta comparison, semantic claim guard, Figure 6 numbers, manuscript-ready summary, and provenance hashes.

Run:

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd <this package>
./run_stage04b8_closeout.sh "$P"
./collect_results.sh
```

Expected final gate:

`PASS_STAGE04B8_FINAL_54G_CLOSEOUT`
