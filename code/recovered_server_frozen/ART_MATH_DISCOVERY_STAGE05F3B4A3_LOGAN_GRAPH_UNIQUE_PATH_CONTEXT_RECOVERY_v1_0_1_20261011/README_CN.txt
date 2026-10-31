Stage05F3B4A3 v1.0.1 — Logan graph unique-path source-context recovery

目的：
B4A2 31/31 boundary anchors 均在 Logan unitig 中 UNIQUE，但 unitig 内 immediate overhang=0。
本阶段只沿 Logan/BCALM2 k=31 assembly graph 的唯一度=1路径恢复缺失 flank。

冻结原则：
- 全部 31 个边界进入；不按 candidate/annotation/结果筛选。
- k=31，邻接 unitig overlap=30 nt。
- 只有 oriented graph degree=1 时才继续。
- 一旦 branch / dead-end / cycle / missing-link，立即停止并记录 censoring。
- 严禁依据 abundance、RT/ART annotation、repeat motif、后续结构结果选择 graph branch。
- 目标仍是 candidate 两侧各 10 kb；达到所需 flank 后停止。
- 本阶段不做 repeat/partner/ART/novelty 判定。

输出：
outputs/F3B4A3_BOUNDARY_GRAPH_EXTENSION_31.tsv
outputs/F3B4A3_RECOVERED_CONTEXT_16.tsv
outputs/contexts_graph_recovered/*.fna
outputs/STAGE05F3B4A3_GATE.json
