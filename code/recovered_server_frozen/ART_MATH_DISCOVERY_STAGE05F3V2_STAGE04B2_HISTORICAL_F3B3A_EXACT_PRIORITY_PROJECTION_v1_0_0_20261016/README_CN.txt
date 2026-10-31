Stage05F3V2 Stage04B2：把 2026-10-11 已冻结的历史 F3B3A compute-priority 规则原样投射到 54-Gb delta 的 DELTA_ONLY combined-family universe。

硬边界：
1) 不改 54-Gb scan / HMM / family clustering；
2) 不使用 known-RT annotation、novelty、structure、truth 做排序；
3) 历史 F3B3A exact script SHA 必须为 d222f638c627566909765012874c7aaa531a1a79b5e67a0d73f726f184a30ed8；
4) 先在历史 6.12-Gb 数据上 byte-identical replay 原 F3B3A outputs；
5) 再从历史 R3 domtblout/R4 map 重建 formal_model_n / qualifying_model_n 语义；
6) 两个审计都 PASS 才允许在 24,810 个 DELTA_ONLY combined families 上执行 exact 4-axis x 8 + deterministic backfill 到 32。

运行：
  conda activate drhu
  bash run_stage04b2_exact_priority_projection.sh /mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1 2>&1 | tee logs/run_stage04b2_console.log

预期最终 gate：
PASS_F3V2_STAGE04B2_HISTORICAL_F3B3A_EXACT_RULE_PROJECTED_TO_DELTA_ONLY_FAMILIES
PASS_STAGE05F3V2_STAGE04B2_HISTORICAL_F3B3A_EXACT_PRIORITY_PROJECTION_COMPLETE
