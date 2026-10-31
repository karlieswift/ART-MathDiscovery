Stage05F3V2 Stage04B3 - DELTA32 Boltz-2 + frozen R5C structure generalization
Date: 2026-10-16

目的：对 Stage04B2 在独立 54.378-Gb delta universe 中按历史 F3B3A 规则冻结得到的 32 个 DELTA_ONLY family representatives，原样复用 2026-10-11 Stage05F3B3C 的 full-3D/R5C protocol。

严格不变：
- shortlist = Stage04B2 exact 32，不允许人工删/换。
- Boltz-2, MSA server, seed=20261007, recycling=3, sampling=200, diffusion_samples=1, --no_kernels。
- 4 个独立 single-GPU workers；不做 DDP/NCCL。
- 不截断长蛋白；primary failure 不允许 msa: empty fallback。
- Foldseek clean-RT R5C: evalue=1.0, max-seqs=100, exact-tmscore=1。
- STRONG: TM>=0.50 AND qcov>=0.50 AND tcov>=0.50。
- SUGGESTIVE: TM>=0.40 AND qcov>=0.30 AND tcov>=0.30。
- technical failure = UNRESOLVED，不得当作 negative。

运行服务器：as-4125
环境：art_structure（Boltz/CUDA）；Foldseek 自动通过 conda run -n art_novelty 调用。

推荐：
  conda activate art_structure
  export BOLTZ_CACHE=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache
  ./run_stage04b3_delta32.sh 32 4 /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1

注意：32 个中存在 >2,000-aa 长序列；不得因为慢或显存压力而截断。失败保留为 technical failure，后续另做资源闭合，不修改 primary protocol。
