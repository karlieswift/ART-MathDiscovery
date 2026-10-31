# ART MathDiscovery Stage05D4-R5B

RT-focused structural reference-panel search for the 43 frozen R5 P1 Boltz-2 full-coordinate structures.

## Frozen reference design

The panel is defined independently of candidate-specific structural results from two external reference arms: curated Swiss-Prot RT-text annotations and the top 10 Swiss-Prot hits per RVT-All 45-profile HMM at full/domain E <= 1e-5. The resulting accessions are intersected with the already frozen AlphaFold/Swiss-Prot Foldseek database, hashed, and converted to a Foldseek structural subdatabase before any P1 search.

## Typical run

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
R4="$P/ART_MATH_DISCOVERY_STAGE05D4R4_PROSTT5_FOLDSEEK_REMOTE_STRUCTURE_v1_0_0_20261007"
R5="$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007"
R5B="$P/ART_MATH_DISCOVERY_STAGE05D4R5B_RT_FOCUSED_STRUCTURAL_PANEL_v1_0_0_20261007"
cd "$R5B"
./run_selftest.sh
RVT=$(find "$P" -type f -name 'RVT-All.hmm' -print | head -1)
echo "$RVT"
./freeze_reference_panel.sh "$P" "$R4" "$RVT" 32
./run_stage05d4r5b.sh "$R5" 32
./collect_stage05d4r5b.sh
```

Do not substitute another HMM file: the script requires SHA256 `74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e`.
