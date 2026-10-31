Stage05F3B3C v1.0.0

目的：对 F3B3B 已冻结的 32-family structure handoff 做 Boltz-2 full-3D，并使用既有 Stage05D4-R5C clean-RT Foldseek DB 与完全相同的冻结规则分类。

固定 Boltz protocol：
- model=boltz2
- MSA server assisted
- seed=20261007
- recycling_steps=3
- sampling_steps=200
- diffusion_samples=1
- --no_kernels
- output_format=mmcif
- 独立单 GPU workers；不使用 DDP/NCCL

固定 R5C：
STRONG: alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50
SUGGESTIVE: alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30
WEAK: 有 hit 但未满足上两档
NO_HIT: 无 Foldseek row

规则：
- 不截断超长蛋白。
- MSA/OOM 等失败记 TECHNICAL_FAILURE；不静默切换 msa: empty。
- known-system matches 保留，不能因 known annotation 删除。
- STRONG/SUGGESTIVE 只表示 frozen clean-RT structural compatibility，不等于 biochemical RT/ART identity。
- WEAK/NO_HIT 不等于 non-RT，也不构成 novelty evidence。
