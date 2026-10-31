F3B2 R2: serialized query-ID -> original MMseqs DB-key mapping hotfix
===================================================================

背景
----
F3B1 的 formal candidate IDs 来自真正送入 hmmscan 的
F3B1_ALL_ORFS_NUMERIC.faa 的 FASTA query headers。

R1 错误地把这些 query IDs 直接当成 original F3B1_orfs_aa MMseqs DB keys
交给 createsubdb，因此只恢复了 10,232 / 10,464。

R2 修复
------
1. 仍冻结并使用原 10,464 formal candidate membership；
2. 顺序扫描 exact hmmscan input FASTA 与 original ORFDB.index；
3. 用 serialized record ordinal 将 hmmscan query ID 映射回 original DB key；
4. 对全部 10,464 条逐条验证 HMM-input peptide == original ORFDB peptide；
5. 直接从 MMseqs sequence/header DB 读取 exact peptide 和 source header；
6. 只有全部 10,464 映射/序列验证 PASS 后才进入 exact-AA/family compression。

不改变
------
- F3B1 formal membership: 10,464
- formal/exploratory HMM thresholds
- RVT-All HMM
- family min-seq-id 0.30
- coverage 0.70
- cov-mode 0
- cluster-mode 0
- sensitivity 7.5
- 不使用 known RT annotation / novelty information

运行
----
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"
tar -xzf ART_MATH_DISCOVERY_STAGE05F3B2_R2_SERIALIZED_QUERYID_KEYMAP_HOTFIX_v1_0_2_20261011.tar.gz
cd ART_MATH_DISCOVERY_STAGE05F3B2_R2_SERIALIZED_QUERYID_KEYMAP_HOTFIX_v1_0_2_20261011
./run_stage05f3b2_r2.sh 32 2>&1 | tee logs/run_stage05f3b2_r2_console.log
