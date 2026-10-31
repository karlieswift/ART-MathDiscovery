Stage05F3B4C：全 10,598 formal RT locus repeat-architecture census + downstream partner candidate scan

目的
1. 对 F3B2-R4 冻结的全部 10,598 formal RVT-supported loci / 4,912 operational families，使用与 F3B4B1 字节完全一致的 repeat detector 做穷尽 census。
2. 不等 repeat 结果再定 partner 规则：同一个 implementation lock 中预先冻结 partner search。
3. 对每个 UNIQUE_EXACT RT locus 的 candidate-oriented downstream DNA 扫描最多 10 kb；code 11；ATG/GTG/TTG start；TAA/TAG/TGA stop；complete ORF >=30 aa；两条 strand 都枚举；物理最近 ORF 作为 primary partner candidate。
4. 最近 downstream ORF protein 用 MMseqs2 做 operational clustering：30% identity / 70% coverage / cov-mode 0 / cluster-mode 0 / sensitivity 7.5。
5. repeat-positive / strict / broad loci 中，若同一 nearest-partner cluster 跨 >=2 独立 accessions 复现，记录为 recurrent nearest-downstream partner cluster。

重要边界
- nearest downstream ORF 是 partner candidate，不等于 confirmed dedicated partner gene。
- repeat geometry + recurrent partner cluster 是 ART-like architecture evidence，不等于 biochemical RT activity / interaction / biological function / novelty。
- no repeat / no partner 只有在 context 足够时才可作局部 observation；不得把 contig truncation 写成 biological absence。
- 不修改 F3B4B1 detector，不修改 F3B2 family membership，不修改 R5C thresholds。

运行
  conda activate art_novelty
  ./run_stage05f3b4c.sh 32 2>&1 | tee logs/run_stage05f3b4c_console.log

结果打包
  ./collect_stage05f3b4c.sh

重点输出
- outputs/STAGE05F3B4C1_GATE.json
- outputs/F3B4C1_GLOBAL_FORMAL_RT_FAMILY_REPEAT_SUMMARY_4912.tsv
- outputs/F3B4C1_ALL_DETECTED_REPEAT_ARRAYS.tsv
- outputs/STAGE05F3B4C2_GATE.json
- outputs/F3B4C2_RT_FAMILY_REPEAT_PARTNER_ARCHITECTURE_SUMMARY_4912.tsv
- outputs/F3B4C2_ARRAY_ASSOCIATED_RT_PARTNER_CANDIDATE_FAMILIES.tsv
