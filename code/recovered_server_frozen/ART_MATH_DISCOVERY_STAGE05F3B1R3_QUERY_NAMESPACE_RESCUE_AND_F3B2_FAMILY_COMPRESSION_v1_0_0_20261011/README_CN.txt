F3B1R3 + F3B2 integrated rescue

原因：F3B1 HMM 输入 FASTA 的 query-name token 并不唯一；R2 在 formal query name=54 处检测到重复。
因此旧 gate 中 formal_candidate_n=10464 只能解释为 namespace-collapsed 的 unique query-name 数，不能直接当成 10,464 个唯一 peptide instance。

R3 修复：
1. 保留冻结的 10,464 formal query-name 集合不变。
2. 扫描原始 65,315,646-record HMM 输入 FASTA，仅提取这些 query name 对应的全部 record。
3. 用全局 FASTA ordinal 赋予每条 peptide 唯一 instance ID。
4. 逐条与原始 MMseqs ORF DB key/sequence 做强校验。
5. 仅对这批 namespace-expanded peptide instances 重跑完全相同的 45-profile RVT-All HMM、完全相同的 1e-4/1e-2 阈值。
6. 要求 corrected formal instances 覆盖原先全部 10,464 formal query names，否则 FAIL。
7. 用 corrected formal peptide instances 做 exact-AA dedupe 和已冻结的 30% identity / 70% coverage operational family compression。

不重跑全 65M HMM；不改变 source universe、RVT-All、阈值、annotation policy 或 family 参数。
