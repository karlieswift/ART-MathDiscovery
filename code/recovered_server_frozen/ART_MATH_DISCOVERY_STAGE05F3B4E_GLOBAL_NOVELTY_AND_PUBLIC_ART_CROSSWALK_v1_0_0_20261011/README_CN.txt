Stage05F3B4E — Global novelty + public ART crosswalk
2026-10-11

目的：尽快补 F3B4D 后剩余的两个计算缺口。
G1) 5 个 RT-fold-supported ART-like family 做 dated public-global sequence/locus audit：
    RT / partner -> NCBI nr remote BLASTP
    frozen architecture core -> NCBI nt remote MEGABLAST
    蛋白 exact/high-near-full 规则继承 F3B3B：99.999/99/99 与 90/80/80。
G2) 对作者公开可精确定位的 8 个 ART source genomes 做直接 crosswalk：
    SA1 MW218148.1 (ART_21)
    MarsHill MW248466.1
    Madawaska MW349129.1 (ART_19)
    LY01 OR836606.1 (ART_20)
    S6 LC680885.1
    PALS_2 MN091626.1
    UFV_DC4 MZ779063.1
    PB50 OR770614.1 (ART_17)

重要边界：
- NCBI nr/nt 是 dated public audit，不等于所有全球 metagenomic/proprietary corpus。
- Anthropic 报告 95 个 90%-identity ART RT clusters，但截至本 freeze 没有可用的完整 machine-readable ART95 sequence catalogue。
- 因此本阶段可以闭合“public NCBI audit + exact public ART8 crosswalk”，不能凭此授权“Anthropic-unreported ART”绝对表述。

运行：
1. conda activate art_novelty
2. ./tests/run_selftest.sh
3. ./run_stage05f3b4e_prepare.sh
4. 两个终端并行：
   ./run_stage05f3b4e_global_remote.sh
   ./run_stage05f3b4e_public_art8.sh
5. ./finalize_stage05f3b4e.sh
6. ./collect_stage05f3b4e.sh

依赖缺失时：
conda install -n art_novelty -y -c conda-forge -c bioconda blast ncbi-datasets-cli
