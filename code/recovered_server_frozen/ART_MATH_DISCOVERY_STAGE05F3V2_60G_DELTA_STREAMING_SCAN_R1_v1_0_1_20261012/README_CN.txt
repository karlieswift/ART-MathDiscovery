Stage05F3V2 60G DELTA STREAMING SCAN R1 v1.0.1 (2026-10-12)

R1 是 acquisition-only hotfix，修复 v1.0.0 的固定 K=5000 容量错误。
已冻结的科学规则完全不变：
- 原 6.121620524-Gb F3 universe 不改、不重算；post-discovery external expansion only。
- source-blind strata 与配额不变：A viral 13.5 Gb；B environmental 18 Gb；C host-associated 13.5 Gb；D extreme 9 Gb。
- source selection 仍只允许 SRA accession/assay/organism/taxonomy + Logan technical file size + realized nucleotide bp；禁止 RT annotation/RVT score/candidate similarity/ACNV/recurrence。
- six-frame maximal stop-to-stop；genetic code 11；no start requirement；min45aa。
- RVT-All SHA256 74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e。
- formal full E<=1e-4 AND domain i-E<=1e-4；exploratory 1e-2/1e-2；family 30% identity/70% coverage/cov-mode0/cluster-mode0/s=7.5。

R1 acquisition changes only:
1. 删除 fixed top-K ceiling。01 对全部 source-blind eligible Logan candidates 使用与 v1.0.0 完全相同的 SHA256(seed|stratum|accession) 64-bit 排序，并用 external sort 冻结完整 deterministic order。
2. 如 v1.0.0 candidate pool 存在，强制验证其每个 stratum 都是 R1 order 的 exact prefix；有任何 drift 即 FATAL。
3. 02 自动导入 v1.0.0 的 work/download/selected.tmp.tsv（或完整 selected manifest），逐条验证已有 .zst SHA256；同一文件系统优先 hardlink 零拷贝复用，跨文件系统退回 symlink；不重新下载。
4. acquisition ledger 记录 SELECTED / LOGAN_OBJECT_NOT_FOUND / HEAD_FAILED_TRANSIENT / DOWNLOAD_FAILED / DECOMPRESS_FAILED / EMPTY_OR_INVALID_FASTA / COMPRESSED_SIZE_OUT_OF_RANGE / REALIZED_BP_OUT_OF_RANGE / ALREADY_IN_OLD_UNIVERSE / DUPLICATE_ACCESSION / ACCEPTED / ACCEPTED_IMPORTED_V100。
5. S3 transport/network failure 与 404 分离。连续 transport failures 达阈值后 fail-closed；以后只需原命令重跑即可继续，不会把网络错误伪装成 source absence。
6. 02 不再因为候选窗口不足而 fatal；它沿完整 deterministic order 消耗候选，直到 A/B/C/D 各自达到 frozen bp quota。只有“完整 eligible universe 真正耗尽仍不足”才报 source-universe limitation。

默认从以下旧目录续跑：
/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_v1_0_0_20261011
可用 F3V2_RESUME_FROM 覆盖。

推荐运行：
  conda activate art_novelty
  P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
  cd "$P"
  sha256sum ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_R1_v1_0_1_20261012.tar.gz
  tar -xzf ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_R1_v1_0_1_20261012.tar.gz
  cd ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_R1_v1_0_1_20261012
  ./tests/run_selftest.sh
  ./01_select_logan_delta_candidates.sh 2>&1 | tee logs/run_01_r1_console.log
  ./02_download_and_freeze_54g_delta.sh 2>&1 | tee logs/run_02_r1_console.log

02 完成且 gate PASS 后，再执行原冻结 scan：
  ./03_run_54g_delta_scan.sh 12 4 2>&1 | tee logs/run_03_scan_console.log
  ./collect_stage05f3v2.sh

重要：不要再运行旧 v1.0.0 的 02_download_and_freeze_54g_delta.sh。旧目录只作为 SHA-verified resume source。
