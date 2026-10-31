Stage05F3B3B known-system audit + 32-family structure handoff

作用：
1. 严格验证 F3B3A frozen shortlist 仍为 32 个 family，候选与 family membership 不改变。
2. 先冻结 F3B3B implementation，再做任何 known-system sequence-similarity 输出。
3. 已知系统审计分两层：
   - ALWAYS: F3B3A 中已经冻结的 best RVT-All profile label（descriptive known-profile annotation）。
   - OPTIONAL ORTHOGONAL: 从项目中既有 pre-F3B frozen reference FASTA 自动建立 known-system protein panel；若无法解析，不阻断 32-family structure handoff，但明确记录 sequence-reference gap。
4. known-system matches 一律保留，不删除 shortlist family。
5. 32 条 structure handoff FASTA 必须与 F3B3A shortlist FASTA byte-exact identical。
6. 本阶段不作 novelty 判定，不改任何 F3B1/F3B2 scientific threshold。

可选环境变量：
F3B3B_REFERENCE_FASTAS=/abs/a.faa:/abs/b.faa
若设置，则只使用这些显式 reference FASTA；否则按 frozen 路径规则扫描既有 ART_MathDiscovery reference panels。
