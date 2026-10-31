Stage05F3V2 60G DELTA STREAMING SCAN v1.0.0

目的：保留原 6.121620524-Gb F3 universe 完全不变，只新增约 54 Gb Logan contig DNA，形成约 60-Gb combined exposure。
这是 post-discovery external expansion，不回写原 F3 primary。

Source selection：只使用 SRA accession / assay / organism / taxonomy 与 Logan technical size / realized nucleotide bp；禁止 RT annotation、RVT score、candidate similarity、ACNV/recurrence 等参与选源。
自动 selector 从 Logan sra_taxid.csv.zst 流式生成四 strata deterministic candidate pool；必须先排除原 F3 至少 131 source accessions。

Delta quotas：A viral-enriched 13.5 Gb；B environmental metagenome 18 Gb；C host-associated metagenome 13.5 Gb；D extreme/underrepresented 9 Gb。总目标 54 Gb；最终 gate 接受 50-60 Gb delta。
Primary sequence representation：Logan contigs。Unitigs 不进入 primary。

Frozen scan：six-frame maximal stop-to-stop；genetic code 11；no start requirement；min45aa；RVT-All SHA256 74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e；formal full E<=1e-4 AND domain i-E<=1e-4；exploratory 1e-2/1e-2；family compression 30% identity/70% coverage/cov-mode0/cluster-mode0/s=7.5。

Minimal run:
  conda activate art_novelty
  ./00_install_or_check_tools.sh
  ./tests/run_selftest.sh
  ./01_select_logan_delta_candidates.sh 5000
  ./02_download_and_freeze_54g_delta.sh
  ./03_run_54g_delta_scan.sh 12 4
  ./collect_stage05f3v2.sh

If tools are missing:
  conda install -n art_novelty -y -c conda-forge -c bioconda awscli zstd hmmer mmseqs2

CPU/IO workload only; no GPU is required.
