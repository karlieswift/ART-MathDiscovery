# R2 v1.0.2 hotfix

在任何 B1 baseline score 产生前修复两项实现语义：
1. Stage02B wide proposal 表不能因为列名同时含 ACNV/recurrence 就把每一行判成 dual-arm；现在 long-format 用显式 family，wide-format 只按冻结 q<=0.005 的 family-specific q 列分类。
2. six-frame / Prodigal HMM hits 在排名前按“同方法+同链+坐标实际重叠(gap=0)”truth-blind collapse 成 RT loci，避免同一基因座多个 peptide hit 人为膨胀 candidate burden。

科学阈值、positive/control panel、RVT HMM、10-kb primary neighborhood、truth chronology 均不变。

# Stage05F4A Formal Baseline Benchmark R1 v1.0.1

这是 v1.0.0 B0 HOLD 后的 **result-blind implementation hotfix**。R1 只修输入解析/历史 artifact 兼容和 runtime 绑定，不改变 benchmark scientific contract。

## R1 修复内容

1. 对 LY01/S6/PALS_2/UFV_DC4/PB50 使用 Stage05A `inputs/genomes/*.frozen.fasta` exact-source precedence，停止从多个副本中猜测。
2. mechanistic controls 绑定到 Stage04C1 实际评分使用的 R3-corrected source：Sen2=`CP001363.1`，DRT2=`AP006725.1`，EC78/Eco7=`CP026202.1`；历史 reference accession 仍单独保留。
3. 支持 Stage05B/Stage04C1 的 `C1_SELECTED_PROPOSALS.tsv` / `C1_UNION_REGIONS.tsv`，并在 B0 truth-blind 地 canonicalize 为统一 proposal/union tables。
4. Prodigal 不再要求当前 conda PATH；B0 自动发现现有 `/home/hgq/conda-envs/*/bin/prodigal`，并冻结 exact path + SHA256。B1 使用冻结的 exact binary path。
5. HMMER 同样使用 B0 冻结的 exact path。
6. RVT baseline 同时要求 full-sequence E-value <=1e-4 与 domain i-Evalue <=1e-4；阈值与 v1.0.0 contract 不变，只修实现与 contract 的一致性。
7. B2 增加 frozen `ARM_M_FULL_FROZEN / ACNV_ONLY_FROZEN / RECURRENCE_ONLY_FROZEN`，因此正式 head-to-head 不再只有 protein-first methods。

## 未改变

- primary positive set / development exclusion / L0050 sensitivity role
- q threshold / recurrence engine / ACNV engine / region merge
- RT neighborhood primary radius ±10 kb；±5/±20 kb 仍仅 sensitivity
- RVT-All HMM SHA256
- truth-blind score-freeze chronology
- no retuning

## 运行顺序

```bash
conda activate art_novelty
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"
sha256sum -c ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_BENCHMARK_R1_v1_0_1_20261013.tar.gz.sha256
tar -xzf ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_BENCHMARK_R1_v1_0_1_20261013.tar.gz
cd ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_BENCHMARK_R1_v1_0_1_20261013
./tests/run_selftest.sh
./00_preflight_and_freeze.sh "$P" 2>&1 | tee logs/run_00_console.log
```

B0 必须得到：
`PASS_STAGE05F4A_B0_R1_INPUTS_BOARDS_TOOLCHAIN_FROZEN_SCORE_READY`

然后才运行：

```bash
./01_run_score_freeze.sh 16 2>&1 | tee logs/run_01_console.log
```

B1 必须得到：
`PASS_STAGE05F4A_B1_R1_PROTEIN_FIRST_BASELINE_SCORES_FROZEN_NO_TRUTH`

B1 结束前不要打开/编辑 truth seal。B2 truth resolution/evaluation 后续单独进行。
