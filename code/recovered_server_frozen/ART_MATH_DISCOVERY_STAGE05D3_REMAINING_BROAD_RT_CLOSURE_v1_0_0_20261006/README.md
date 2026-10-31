# ART MathDiscovery Stage05D3 remaining broad-RT closure

Purpose: audit the three Stage05D1 formal broad-RVT targets not included in Stage05D2: MT708547.1, NC_031049.1, and PV251870.1.

## Windows PC: retrieve official GenBank records
Run `windows_download_official_genbank.ps1` from PowerShell in a clean directory, then copy the three `.gb` files to this package under `inputs/official_pc/` on the analysis host.

## Analysis host
```bash
conda activate art_novelty
D1=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D1_NOVELTY_FIRST_MINING_v1_0_0_20261006
./run_selftest.sh "$D1"
./run_stage05d3.sh "$D1"
./collect_stage05d3.sh
```

Expected final gate: `PASS_REMAINING_STAGE05D1_BROAD_RT_HITS_AUDITED_WITH_NOVELTY_HOLD`.
