Stage05F3B3E — serial-exclusive resource-only technical closure

目标：对 F3B3D exact-protocol retry 后仍 unresolved 的全部且仅限 7 条 Boltz CUDA-OOM 候选，做最后一次 resource-only technical closure。

唯一执行层变化：
- 之前：4 个独立单 GPU worker 并发。
- 本阶段：选一张无 compute process 的 GPU；7 条严格串行；每条使用新的 Boltz subprocess。

科学设置完全不变：
- boltz2
- MSA server-assisted
- seed 20261007
- recycling 3
- sampling 200
- diffusion samples 1
- --no_kernels
- --devices 1
- GPU
- mmCIF
- 不截断
- 不允许 msa: empty fallback
- R5C Foldseek 与 STRONG/SUGGESTIVE threshold 不变

本阶段是 post-hoc resource-only technical closure，不替代 F3B3C primary 或 F3B3D exact-protocol retry。
若仍失败，只能标记 resource-limited unresolved，不能解释为结构阴性。

运行：
  conda activate art_structure
  export BOLTZ_CACHE=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache
  ./run_stage05f3b3e.sh 32 auto

第二参数可显式指定物理 GPU，例如：
  ./run_stage05f3b3e.sh 32 2
