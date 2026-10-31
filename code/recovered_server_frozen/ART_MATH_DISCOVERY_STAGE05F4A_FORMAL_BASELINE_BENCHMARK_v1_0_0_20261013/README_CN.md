# Stage05F4A Formal Baseline Benchmark v1.0.0

## 目标

这是**论文级 formal baseline comparison**，不是为了让我们的结果“赢”而调参。所有 baseline 在 truth 揭盲前冻结，primary comparison 结果无论方向如何都保留。

### Primary external positives（8）
Madawaska、Machias、SA1、LY01、S6、PALS_2、UFV_DC4、PB50。

- MarsHill = development only，不进入 primary inferential comparison。
- L0050 = locked retrospective sensitivity，不进入 primary external count。
- Sen2 / Ec78(Eco7) / DRT2 = mechanistic non-array RT controls。

## Primary comparison

1. `ARM_M_FULL_FROZEN`：已有 frozen ACNV-family + recurrence whole-genome board；不重算。
2. `PRODIGAL_RVT_PLUS_FROZEN_REC`：传统 protein-first：Prodigal meta -> myRT RVT-All -> 只在 RT ±10 kb 邻域读取**同一 frozen recurrence proposal board**。
3. `SIXFRAME_RVT_PLUS_FROZEN_REC`：更强的 annotation-independent protein-first upper bound：six-frame stop-to-stop >=45 aa -> RVT-All -> 同一 frozen recurrence board。

这样可以公平隔离：优势究竟来自 raw-DNA-first，还是只是 repeat detector 更强。

## Secondary/decomposition

- ACNV-only frozen
- recurrence-only frozen
- RVT-only（Prodigal / six-frame）
- TRF repeat-only 与 TRF+RVT 为 optional secondary；TRF 缺失不会阻止 primary benchmark。

TRF 的限制必须在论文中明确：TRF 定义的是 adjacent approximate tandem repeats，而 ART 的短 repeat core 可被长、非同源 spacer 隔开，因此不能把 TRF 当成唯一 conventional baseline。

## 冻结规则

- RVT-All HMM 必须匹配项目已冻结 SHA256：`74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e`。
- myRT threshold：`E <= 1e-4`, `domE <= 1e-4`。
- Prodigal primary：`-p meta`。
- six-frame ORF：stop-to-stop，最短 45 aa。
- RT neighborhood primary：±10 kb；±5 kb / ±20 kb 仅 sensitivity。
- recurrence：不重算；读取已有 frozen proposal board。
- `01_run_score_freeze.sh` 在代码层面不能读取 truth seal。
- score freeze 之后才执行 `02_postfreeze_unblind.sh`。
- `retuning_authorized=false`。

## 第一步：只做 B0 preflight

```bash
conda activate art_novelty

P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
cd "$P"

sha256sum -c ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_BENCHMARK_v1_0_0_20261013.tar.gz.sha256

tar -xzf ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_BENCHMARK_v1_0_0_20261013.tar.gz

cd ART_MATH_DISCOVERY_STAGE05F4A_FORMAL_BASELINE_BENCHMARK_v1_0_0_20261013

./tests/run_selftest.sh
./00_preflight_and_freeze.sh "$P"

cat frozen/BASELINE_B0_GATE.json
column -t -s $'\t' frozen/BASELINE_SEQUENCE_MANIFEST.tsv
column -t -s $'\t' frozen/BASELINE_BOARD_INVENTORY.tsv
cat frozen/BASELINE_TOOLCHAIN.json
```

**先不要直接跑 01，除非 B0 gate 是 PASS。**

B0 会 fail closed，而不是猜路径。它会：
- 找 13 个 frozen genome/contig FASTA；
- 若同 accession 有多个 copy，要求 normalized sequence SHA256 一致；
- 按 frozen SHA 找 RVT-All.hmm；
- 检查 Prodigal/HMMER；
- inventory 现有 frozen proposal/union boards；
- 搜索可能的 post-freeze truth source files，但不打开 truth 做 score。

## B0 PASS 后

```bash
./01_run_score_freeze.sh 16
cat outputs/score_freeze/BASELINE_B1_GATE.json
```

B1 只产生 protein-first candidate scores；不读取 truth。

## truth 解封

当前 v1.0.0 已直接冻结 Madawaska/Machias/SA1/L0050/PB50/DRT2 的已有坐标；LY01/S6/PALS_2/UFV_DC4 以及 Sen2/Ec78 必须从项目既有 post-freeze truth artifacts 精确恢复。B0 会输出：

`outputs/TRUTH_SOURCE_CANDIDATES.tsv`

只有把这些坐标绑定到既有 artifact SHA 后，才更新一份新的 truth TSV 作为 `02_postfreeze_unblind.sh` 输入；不得看 B1 结果后选择更有利的坐标。

## 最终论文 endpoint

Primary positive set n=8：
- recovered
- truth rank / candidate N
- MRR
- Hit@1 / @3 / @5 / @10
- normalized rank percentile
- candidate burden / Mb
- RT-neighborhood exposed genome fraction

Controls n=3：
- known RT-system locus 是否进入 recurrence-supported RT-first candidate
- 与 frozen Arm M `0/3 dual-arm` 对照

统计以 paired descriptive metrics、exact paired tests 和 bootstrap CI 为主；不对 n=8 做夸张的 asymptotic claim。
