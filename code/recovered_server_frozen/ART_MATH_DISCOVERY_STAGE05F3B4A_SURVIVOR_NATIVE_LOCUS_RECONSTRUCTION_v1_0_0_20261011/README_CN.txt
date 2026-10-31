Stage05F3B4A — 16 structural-survivor native-locus reconstruction

目的：
- 只对 F3B3D 已冻结的 STRONG/SUGGESTIVE 16 个 structural survivors 做 source-contig exact placement。
- 从 F3B2-R4 canonical source mapping + original MMseqs nucleotide DB 恢复 exact parent contig。
- 通过六框 translation exact peptide replay 定位 candidate ORF；不依赖 annotation，不猜坐标。
- 每个 candidate 输出 genomic ±10 kb context 与 candidate-oriented ±10 kb context。
- 不做 repeat/partner/ART/novelty 结论；这里只生成下一阶段的冻结 locus input。

硬边界：
- survivor membership 不变；必须来自 F3B3D final STRONG/SUGGESTIVE。
- 预期 16/16。
- exact peptide placement 必须唯一；0 hit 或 >1 hit 均 HOLD，不猜。
- genetic code = 11；start-codon requirement = none；与 F3B1 translation science 一致。
- flank = 10,000 bp each side；只是 context materialization，不是 biological threshold。

推荐：r750-1 CPU 跑。GPU 不需要。
