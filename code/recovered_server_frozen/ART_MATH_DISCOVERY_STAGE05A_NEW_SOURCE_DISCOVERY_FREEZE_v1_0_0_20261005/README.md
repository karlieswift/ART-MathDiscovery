# Stage05A：新来源噬菌体发现队列

## 目标与范围

本包推进“发现新的 ART-like 系统，筛出湿实验候选”的主线：从 NCBI Datasets 获取此前未纳入本项目的完整 Caudoviricetes DNA 来源，分层选择、下载并冻结。它不是重跑原240个genome的审计，也不把已知RT当作新ART。

**本次新增：可执行的 metadata 获取、固定选择、完整 genome 下载、序列去重、冻结和计算规模计划。**

**本包不执行数学评分、ORF、HMMER、array搜索或wet lab。** 源集合先锁定，按真实规模安排下一次隔离的whole-genome两臂评分；不得根据评分结果改变入组。UG27、原Stage04A、383候选、R1D2和R1E均不读取修改。

## 已实测／未执行

已执行25项离线自测，覆盖metadata双schema、拒绝冲突字段、来源过滤、确定性配额、序列解析、MD5、反向互补、N排除及原240排除清单。排除清单使用真实原240序列构建，共16,574,805 bp。

本端没有执行新的NCBI在线下载，没有得到实际新增1000个genome，也没有产生新增ART候选数量。服务器实际结果决定分母。这里的1000是分层配额上限，不是预报取得数量。

## 输入选择规则

- 完整Caudoviricetes记录，要求明确accession.version和所需谱系。
- 有Bacteria taxid2的宿主谱系；或者宿主缺失、但病毒名称包含phage，后者明确标记。已报告却不能确认细菌的宿主不进入本轮。
- 20–100 kb：最多250；100–200 kb：最多250；200 kb–2 Mb：最多500。
- 每个报告的viral taxid保留一个代表。优先RefSeq，其次固定SHA256顺序；这不等于物种独立或近同源去重。
- 排除原240 accession root和7个登记参考来源；再用原240的真实序列和反向互补SHA256排除重复。没有宣称排除了世界上全部已知ART；系统级新颖性要在后续完成。
- 在每层按固定哈希排序取样。不按RT名称、注释是否完整或候选得分选择。
- 下载后N或其他非ACGT序列、多记录FASTA按预先规则排除，不修碱基、不拼接、不重新补样；长度/版本不符或下载未解决则HOLD。
- 精确／反向互补重复排除不替代近同源、系统簇独立性分析。

新增到项目≠生物学新系统≠发表日期前瞻队列。来源快照记录真实取得时间，releaseDate原样记录。

## 推荐运行环境与命令

只需要Python 3.8+标准库和现成的NCBI datasets CLI。使用此前已成功运行datasets的r750-1 / drhu，不需要改bioh200或UG27环境。

```bash
conda activate drhu
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"
sha256sum -c ART_MATH_DISCOVERY_STAGE05A_NEW_SOURCE_DISCOVERY_FREEZE_v1_0_0_20261005.tar.gz.sha256
tar -xzf ART_MATH_DISCOVERY_STAGE05A_NEW_SOURCE_DISCOVERY_FREEZE_v1_0_0_20261005.tar.gz
cd ART_MATH_DISCOVERY_STAGE05A_NEW_SOURCE_DISCOVERY_FREEZE_v1_0_0_20261005
./run_selftest.sh
./run_stage05a_source_freeze.sh
```

先单独查看metadata选择，可运行：
```bash
python3 scripts/discovery_source.py select
```
metadata与selection第一次生成后锁定；续跑不重新拉取metadata改变样本。如果datasets不存在，先使用已安装工具的机器运行，不会静默改装现有环境。

## 断点与故障

重新运行 `./run_stage05a_source_freeze.sh` 复用已验证ZIP和冻结metadata。下载按每条最多3次执行，并记录日志；失败不被当作无RT或非ART。

不要同时启动两个该包任务。`.RUNNING`记录当前PID；意外断电后确认进程已不存在，才可手工移除本包的锁文件。不得删除原项目或UG27目录。

若一个现存ZIP损坏，程序停止而不静默覆盖：先保留故障文件用于检查。`SOURCE_ACQUISITION_LEDGER.tsv`逐条记录排除或未解决原因。

## 输出

- `outputs/SELECTION_FREEZE.json`：metadata快照、选择清单哈希、实际配额与短缺。
- `outputs/SELECTED_SOURCE_ACCESSIONS.tsv`：固定入选的metadata来源。
- `outputs/FROZEN_NEW_GENOME_MANIFEST.tsv`：可评分的实际DNA集合、原序列和FASTA字节哈希。
- `outputs/SOURCE_ACQUISITION_LEDGER.tsv`：逐来源下载／QC／重复排除状态。
- `outputs/SOURCE_COHORT_FREEZE.json`：只有完整处理完所选来源才放行。
- `outputs/COMPUTE_RESOURCE_PLAN.json`：实际窗口量与累计数组量；不是RAM峰值或CPU时长预测。

成功gate：`PASS_NEW_SOURCE_COHORT_FROZEN_WITH_ACTUAL_DENOMINATOR`。即使配额不足也按真实集合冻结，不宣称1000/1000。

## 回传

无论PASS或HOLD，运行：
```bash
./collect_results.sh
```
默认打包代码、日志、metadata、清单和gate，不重复上传全部genome和ZIP。确实需要所有新序列时可显式运行 `./collect_results.sh --with-genomes`。输出包和.sha256一起回传。

## 来源

NCBI Datasets官方命令与schema：
- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/get-virus-metadata/
- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/summary/virus/genome/datasets_summary_virus_genome_taxon/
- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/virus/
- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/virus/genome/datasets_download_virus_genome_accession/

本包不要求E-utilities；沿用该项目已经成功的datasets取数通道。