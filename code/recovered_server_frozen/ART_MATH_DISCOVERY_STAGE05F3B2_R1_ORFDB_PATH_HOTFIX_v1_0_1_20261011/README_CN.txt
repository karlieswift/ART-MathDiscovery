F3B2 R1 ORFDB PATH HOTFIX
========================

问题：
F3B1 v1.0.1 的 work/ 目录只保留了 F3B1_ALL_ORFS_NUMERIC.faa 的 symlink，
没有 F3B1_orfs_aa MMseqs database，因此 createsubdb 指向 S1/work/F3B1_orfs_aa 会失败。

正确来源：
ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010/work/F3B1_orfs_aa

这是生成 F3B1_ALL_ORFS_NUMERIC.faa 的原始 extractorfs MMseqs DB。
F3B1 v1.0.1 的 64-way sharding只改变 hmmscan 调度，没有重建 ORF universe。

R1 保持不变：
- frozen formal candidates = 10,464
- exact-AA dedupe rule
- family min-seq-id = 0.30
- coverage = 0.70
- cov-mode = 0
- cluster-mode = 0
- sensitivity = 7.5
- 不使用 known RT annotation
- 不做 novelty screen

安装/运行：
  conda activate art_novelty
  P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
  cd "$P"
  tar -xzf ART_MATH_DISCOVERY_STAGE05F3B2_R1_ORFDB_PATH_HOTFIX_v1_0_1_20261011.tar.gz
  cd ART_MATH_DISCOVERY_STAGE05F3B2_R1_ORFDB_PATH_HOTFIX_v1_0_1_20261011
  ./run_stage05f3b2_r1.sh 32 2>&1 | tee logs/run_stage05f3b2_r1_console.log

预期最终 gate：
PASS_F3B2_FORMAL_CANDIDATE_FREEZE_AND_FAMILY_COMPRESSION
PASS_STAGE05F3B2_R1_ORFDB_PATH_HOTFIX_AND_FAMILY_COMPRESSION_COMPLETE
