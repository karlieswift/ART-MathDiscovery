# Stage05F4A Formal Baseline — B2 Closeout R4 v1.0.4

这是 **post-score-freeze / B2-only** 包。不要重新运行 B0/B1。

目的：
1. 从已归档的 R2 B1 score-freeze tar.gz 原样导入 `frozen/` 与 `outputs/score_freeze/`；
2. 在 B1 冻结之后，确定性解析此前 seal 中的 `PROJECT_TRUTH_RESOLVE`；
3. 用预注册的 RT-centered ±10 kb neighborhood 评价 protein-first localization；
4. 补齐 contract 已预注册的 non-array RT control dual-signal endpoint；
5. 保留预注册 8-case primary positive panel，同时分层输出：3-case strict external array subset、4-case project-reconstructed public-source array subset、7-case exact-array combined descriptive subset；PB50 仍是 prelocked RT-proxy challenge。

## Truth resolution rule

- LY01/S6/PALS_2/UFV_DC4：从既有 Stage05C2B `CORE_CANDIDATE_QUALIFYING_REPEAT_DESCRIPTIONS.tsv` 中，对该 public sequence 的全部 `EXACT_CORE_PRIMARY_DESCRIPTION` 行取 `min(start0)..max(end0)`，不读取本次 baseline 的任何 score/rank。
- Sen2/Ec78/DRT2：从既有 Stage04C1B0 `C1B0_CONTROL_TRUTH_RESOLVED.tsv` 复制冻结 RT anchor。
- Madawaska/Machias/SA1 与 PB50：保持 CASES 中预锁 interval 不变。

这些 project-reconstructed exact-repeat truth 与真正独立 external truth 的证据层级不同，因此输出中保留 `truth_independence_tier`，不能把它们写成“完全独立 external annotation”。

## Run

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"
sha256sum -c ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_B2_CLOSEOUT_R4_v1_0_4_20261013.tar.gz.sha256
tar -xzf ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_B2_CLOSEOUT_R4_v1_0_4_20261013.tar.gz
cd ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_B2_CLOSEOUT_R4_v1_0_4_20261013
./tests/run_selftest.sh

B1ARC="$P/ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_BENCHMARK_R2_v1_0_2_20261013/ART_MATH_DISCOVERY_STAGE05F4A_B1_R2_SCORE_FREEZE_20261013_172035.tar.gz"
./00_import_b1_score_freeze.sh "$B1ARC"
./01_build_truth_seal.sh "$P"
```

先检查 truth seal 与 gate；确认后才运行：

```bash
./02_postfreeze_unblind.sh
```
