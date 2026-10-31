Stage05F3B3E-R1：Foldseek-only technical repair
- 只处理 F3B3E 中 5 个 Boltz 已成功但 Foldseek rc=1 的 CIF。
- 不重跑 Boltz；不改结构；不改 R5C 数据库、阈值或 Foldseek科学参数。
- 运行前必须 `conda activate art_novelty`，本包直接调用该环境中的 foldseek，故意绕开此前 `conda run -n art_novelty` handoff。
- 若仍失败，日志会直接显示最后 20 行；不得把失败解释为 non-RT。
