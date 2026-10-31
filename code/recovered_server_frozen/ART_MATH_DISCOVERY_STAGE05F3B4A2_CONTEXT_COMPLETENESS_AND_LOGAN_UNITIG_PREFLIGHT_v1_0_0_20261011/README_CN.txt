Stage05F3B4A2 — context completeness audit + Logan unitig preflight

为什么插入这一阶段：
B4A 16/16 exact placement 虽然 PASS，但 source parent contigs 可能远短于 ±10 kb。
在 source context 不完整时，不能把“未看到 repeat/partner/array”解释为 biological absence。
本阶段先量化 contig-edge censoring，并冻结全部 12 source accessions；不删任何 survivor。

默认 preflight：只做 context audit + 检查 Logan unitig URL/size，不下载大文件。
  ./run_stage05f3b4a2.sh preflight

若 unitig 文件规模可接受，再运行：
  ./run_stage05f3b4a2.sh download_scan 2

第二模式会下载全部 12 source-accession unitig files，并精确扫描缺失 flank 的 parent-contig boundary 31-mers。
这里只判定 boundary anchor 是否 UNIQUE/AMBIGUOUS/ABSENT 和 unitig 内立即可恢复的 overhang；不做 biological filtering、不做 array absence claim。
后续只有在 source context 恢复/分支状态明确后，才进入 B4B array/partner architecture。
