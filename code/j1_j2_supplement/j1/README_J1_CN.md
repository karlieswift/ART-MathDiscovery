# J1：原始 clean-RT 校准的可重复性补充

本模块包含可独立检查的原始 R5C 科学资料、已完成的逐行复核结果，以及从原项目只读提取数据库、制作小型可移植参考库并重新运行原协议的脚本。它不修改稿件、图表、旧阶段数据或旧结果目录。

## 已完成的真实复核

原始资料来自 2026-10-07 冻结的 R5C 结果档案；86 份关键科学文件的身份由 `reference/R5C_ARCHIVE_IDENTITY.json` 固定，其整体 SHA-256 又固定在程序中。程序同时检查原档案自带的校验清单。修改原始资料并同步改写其自带清单，也不能被认定为原档案复现通过。

| 原校准面板 | 原脚本：最高 TM hit | 统一规则：所有合格 hits 的最高 tier |
|---|---|---|
| 32 个 known RT，排除同 accession 自身命中 | 29 Strong、1 Suggestive、2 Weak | **31 Strong、1 Suggestive** |
| 96 个长度匹配 non-RT controls | 43 Weak、53 No-hit | **43 Weak、53 No-hit** |

因此，统一规则下正例 Strong + Suggestive 为 32/32，但 Strong 是 31/32。后续冻结文件中的“32 个全为 Strong”与原始数据不符。原 R5C 汇总的 0.9375 对应原算法的 30/32，已逐行复现。三条变化为 P08548：Weak→Strong；P19593：Suggestive→Strong；P31843：Weak→Suggestive。

校准共 1,126 行原始比对：positive 1,013 行，排除 32 条自身比对后剩 981 行；negative 113 行。53 个没有 `.m8` 命中行的 negative controls 均在原结果数据库索引中有独立键，且相应记录为一个 NUL 字节。程序检查所有分片、索引覆盖、偏移、记录末尾、每条 query 的 target-key 多重集及命中行数，所以这里的 53 个 No-hit 有实际完成证据。缺失文件、空文件和格式无效文件分别处理。

32 条自身比对的原始 TM 值为 1.001–1.007，程序原样记录并排除，不截断为 1、不用于非自身校准。非自身比对出现非有限或超范围分数会导致检查失败。

原控制设计仅排除了自身 accession；没有按序列相似性或家族进行留出。96 个 controls 均与对应 RT 长度完全相同，由原 Swiss-Prot 描述排除规则选择，不能改称经湿实验确认的阴性或独立家族外部验证集。0/96 是这组控制上的观测结果。

`local_results/original_calibration_audit/` 提供 128 行控制结果、全部 hit 行、3 条差异、新旧代表命中及原文件行号、53 条空结果证据、原始哈希核对和完整 JSON 报告。`local_results/target_rule_consistency/` 用同一解析器和分级函数重查了现有 114 个 target memberships 中可读取的 112 份结果、2,490 个 hits。所有结果一致；两例历史技术失败继续保留为未解决。strict 54-Gb 仍为 18 Strong、6 Suggestive、7 Weak、1 No-hit。

## 命令接口

总包的一键入口由上层 README 提供。下面是本模块独立接口，所有新输出必须位于项目物理目录内的 `exports` 或 `outputs`。默认使用 `exports`；已有非空输出目录不会覆盖。输出目录也不能嵌入原始校准资料、收集资料或 portable 数据源中。

```bash
ART_PROJECT="/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1"
ART_KIT="$ART_PROJECT/tools/ART_MATH_DISCOVERY_J1J2_COMPUTATIONAL_SUPPLEMENT_v1_0_0_20261022"
ART_RUN="$ART_PROJECT/exports/J1_$(date -u +%Y%m%d_%H%M%S)"

python "$ART_KIT/j1/j1_reproducibility.py" collect \
  --project-root "$ART_PROJECT" \
  --output "$ART_RUN/collection"

python "$ART_KIT/j1/j1_reproducibility.py" audit-calibration \
  --project-root "$ART_PROJECT" \
  --source "$ART_KIT/j1/reference/r5c_archived_results" \
  --output "$ART_RUN/archived_audit"

python "$ART_KIT/j1/j1_reproducibility.py" plan-replay \
  --project-root "$ART_PROJECT" \
  --collection "$ART_RUN/collection" \
  --source "$ART_KIT/j1/reference/r5c_archived_results" \
  --threads 32 --output "$ART_RUN/plan"

python "$ART_KIT/j1/j1_reproducibility.py" export-db \
  --project-root "$ART_PROJECT" \
  --plan "$ART_RUN/plan/REPLAY_PLAN.json" \
  --output "$ART_RUN/portable"
```

`export-db` 已包含在 portable 库上重新运行原 128-control 协议的步骤。如还要在原服务器的原数据库上单独重跑：

```bash
python "$ART_KIT/j1/j1_reproducibility.py" run-replay \
  --project-root "$ART_PROJECT" \
  --plan "$ART_RUN/plan/REPLAY_PLAN.json" \
  --output "$ART_RUN/original_database_replay"
```

`collect` 自动定位的阶段名必须完全等于：

```
ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007
```

没有匹配阶段或有多个候选会返回 `NOT_COMPLETE`，不会选择日期最新的目录。确实需要指定位置时，使用 `collect --r5c-dir /完整/原阶段/路径`。

## 收集、完整性与数据库大小

`collect` 保留阶段内小型科学文件的相对路径，包括零字节 `.m8`、原脚本、配置、日志、原名单、序列、原结果数据库小分片和索引。默认单文件复制上限为 64 MiB；超过上限会进入待补清单，不会伪装成已完成。

对 `references/clean_rt_afsp*` 和 `references/matched_nonrt_afsp*` 的每个实际组件，包括符号链接目标，程序完整计算 SHA-256 并记录大小、解析后的真实路径和时间信息。相同真实目标会复用同一次哈希，避免两组子库对共享 AlphaFold 大文件反复扫描。**原大数据库载荷不会复制进回传目录。** 初次全量哈希可能主要消耗磁盘读取时间；日志会显示当前正在读取的组件。

执行前后还检查完整组件集合、真实目标、大小、mtime、ctime、inode 和 device。新增 `.idx` 等组件、替换符号链接、改写文件后恢复 mtime 等情况都会阻止继续。这里采用“一次完整哈希加执行前后文件状态检查”，不会把它描述为每个阶段都重新读取所有大型载荷。

历史四份清单只记录六个 index/type 组件。六项一致只能说明这六项一致，不能单独证明历史完整数据库身份。`DB_FULL_COMPONENT_SNAPSHOT.json` 与 `NONRT_DB_FULL_COMPONENT_SNAPSHOT.json` 对“当前完整载荷快照”和“历史完整载荷身份”作了明确区分。

## 小型 portable 库

`export-db` 使用已安装的原版本 Foldseek，分别对 AA、3Di、C-alpha 和 header 四种组件调用官方 `base:createsubdb --subdb-mode 0 --id-mode 0`。每个组件在不同的新目录中生成，随后只复制其普通核心载荷、index 和 dbtype 文件。程序不会跟随工具附带的 ancillary 符号链接复制整个 Swiss-Prot/AlphaFold 母库，也不会手写或改造 MMseqs 二进制格式。

完成条件包括：

1. 各组件键集合与原 32/96 清单完全一致，输出无外部符号链接。
2. 官方 `convert2fasta` 导出的源库和 portable 库名称、完整氨基酸序列及原名单长度逐条一致。
3. 完整记录 portable 组件哈希。
4. 仅使用 portable 库重跑原 `search` + `convertalis` 参数。全部 128 个控制的 eligible hit 多重集逐数值一致，才返回 `PASS_PORTABLE_DB_AND_EXACT_REPLAY`。

只复现 tier、但命中值存在差异时会返回待审核状态；不会称为逐值重现。`.m8` 文件行顺序不同而规范化命中多重集完全一致，可视为语义一致，同时报告文件本身是否逐字节相同。

回传应包含完整 `portable` 导出目录中的报告、`reference_portable`、两组源/portable FASTA、128 行身份表、portable replay 的原始结果及日志。`work` 是本次运行的派生中间文件，可由上层压缩打包器排除；原阶段不会被清理或覆盖。

把完整导出目录移到另一台机器后，可用下面接口重新验证与搜索，原服务器绝对路径无需存在：

```bash
python "$ART_KIT/j1/j1_reproducibility.py" replay-portable \
  --project-root "$ART_PROJECT" \
  --portable-root "$ART_PROJECT/exports/received_portable" \
  --source "$ART_KIT/j1/reference/r5c_archived_results" \
  --threads 32 --output "$ART_PROJECT/exports/portable_independent_replay"
```

## 版本、规则与仍待完成的部分

原日志记录 Foldseek `10.941cd33`；官方相同 release 的 `version` 命令返回完整提交号 `941cd33ff0771cd2e3f144e3293e22a2b87e9fda`。程序接受这两种同源标识，记录实际可执行文件哈希；不同版本会暂停，不会静默替换。

原 128-control 搜索是数据库对数据库的 `search` 后接 `convertalis`，采用 E-value 1.0、max-seqs 100、exact-tmscore 1 和原 11 列输出。它与 target 的 `easy-search` 包装形式不同，因此本模块保留原控制协议。分级函数不增加额外 E-value/bit-score 筛选。Strong 门槛为 TM≥0.50 且双向覆盖均≥0.50；Suggestive 门槛为 TM≥0.40 且双向覆盖均≥0.30。每条 query 取所有符合条件命中的最高 tier，并在该 tier 内选择代表命中。

随包 `tests/TEST_REPORT.json` 给出实际测试状态。包含真实官方 Foldseek 的小型结构夹具集成：创建母库、mode-1 子库、四组件 mode-0 portable 导出、序列比对、原/portable 搜索逐字节一致。这是软件接口验证，**不是**原 128-control 数据库搜索已经运行的证据。

截至随包本地审计完成时，原 128 个控制的 raw-hit 重聚合与 114 target 规则一致性检查已完成；原服务器完整数据库的 fingerprint、portable 物化与原协议新搜索仍须由上述命令产生实际结果。没有这些运行产物时，不能把 J1 整体标为已关闭。

## 官方实现依据

- Foldseek 官方结构子库命令：[FoldseekBase.cpp](https://github.com/steineggerlab/foldseek/blob/master/src/FoldseekBase.cpp)。
- 官方 MMseqs2 参数明确区分 mode 0 复制数据、mode 1 软链接载荷：[Parameters.cpp](https://github.com/soedinglab/MMseqs2/blob/master/src/commons/Parameters.cpp)。
- 官方实现仍可能为 ancillary 文件创建链接，因此本模块仅装配普通核心组件：[createsubdb.cpp](https://github.com/soedinglab/MMseqs2/blob/master/src/util/createsubdb.cpp)。
- 最终可用性以原 release 二进制的实际 help、真实结构夹具运行及生产数据执行日志共同核实，不依据最新文档推定旧版行为。
