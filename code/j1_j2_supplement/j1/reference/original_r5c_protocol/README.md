# ART MathDiscovery Stage05D4-R5C

Fast structural calibration after R5B.

## Run

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
R4="$P/ART_MATH_DISCOVERY_STAGE05D4R4_PROSTT5_FOLDSEEK_REMOTE_STRUCTURE_v1_0_0_20261007"
R5="$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007"
R5B="$P/ART_MATH_DISCOVERY_STAGE05D4R5B_RT_FOCUSED_STRUCTURAL_PANEL_v1_0_0_20261007"
R5C="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007"
cd "$R5C"
./run_selftest.sh
./freeze_stage05d4r5c.sh "$R4" "$R5B" 32
./run_stage05d4r5c.sh "$R5" 32
cat outputs/STAGE05D4R5C_SUMMARY.json
cat outputs/GATE_STAGE05D4R5C.json
./collect_stage05d4r5c.sh
```

This stage is calibration/interpretation only. No DNA score, DNA rank, candidate universe, or P1 selection is changed.
