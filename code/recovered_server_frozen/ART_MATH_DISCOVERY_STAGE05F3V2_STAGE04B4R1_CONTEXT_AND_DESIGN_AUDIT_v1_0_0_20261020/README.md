# Stage05F3V2-Stage04B4R1 — 12-representative context-completeness + design audit

Run the lightweight audit first (no unitig download):

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"
tar -xzf ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4R1_CONTEXT_AND_DESIGN_AUDIT_v1_0_0_20261020.tar.gz
cd ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4R1_CONTEXT_AND_DESIGN_AUDIT_v1_0_0_20261020
./run_stage04b4r1.sh "$P" head
```

`head` performs HTTP HEAD requests for the 12 Logan unitig files and writes a download-size plan, but does not download them.

Inspect:

```bash
column -t -s $'\t' outputs/DELTA12_REP_CONTEXT_COMPLETENESS.tsv
cat outputs/SIX_VS_54G_DESIGN_AUDIT.json
column -t -s $'\t' outputs/DELTA12_UNITIG_DOWNLOAD_PREFLIGHT.tsv
```

If the 12 unitig downloads are acceptable, download them:

```bash
bash outputs/DOWNLOAD_DELTA12_UNITIGS.sh unitigs
```

Then run exact terminal-anchor search:

```bash
./run_stage04b4r1.sh "$P" local "$PWD/unitigs"
column -t -s $'\t' outputs/DELTA12_UNITIG_ANCHOR_MATCHES.tsv
```

Do not choose a graph branch manually after seeing results.
