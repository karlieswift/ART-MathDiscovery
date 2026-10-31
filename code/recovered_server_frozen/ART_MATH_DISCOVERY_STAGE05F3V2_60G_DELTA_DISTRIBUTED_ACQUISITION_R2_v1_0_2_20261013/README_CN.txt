F3V2 54G Distributed Acquisition R2 v1.0.2 (2026-10-13)

目的
====
仅把已经冻结的 A/B/C/D 四个 source strata 的 acquisition 分布到多台主机并行执行。
科学设计不变：候选 universe、hash order、每层 bp quota、compressed-size 规则、realized-bp 规则、旧 universe 排除、后续 six-frame/HMM/MMseqs/repeat/structure/novelty 阈值均不变。

为什么安全
========
R1 本来就是 A/B/C/D 四层独立 quota，并且每层内部按照固定 deterministic_rank 顺序扫描，达到本层 bp target 即停止。R2 只让四层同时运行；同一层仍严格 rank-order，因此最终选样等价于串行 R1 的 A/B/C/D 选择（运行顺序不同，不改变每层成员）。

绝对不要
========
1. 不要让旧 R1 的 02_download_and_freeze_54g_delta.sh 与 R2 同时写状态。切换前先 Ctrl-C 停旧 02；已接受文件不会丢。
2. 不要在两个主机上同时跑同一 stratum；R2 有 stratum lock 并会拒绝。
3. 不要把单一 stratum 随意切 rank ranges 后“谁先下载完就接受谁”；那会改变 stop-at-quota 的 deterministic selection。

推荐三机分工
============
as-4125 : A
r750-1  : B
r750-2  : D  (先跑 D，因为候选池只有 18,881，可尽早判断 9 Gb quota 是否可满足)
D 完成后 r750-2 再跑 C。
也可在确认网络/磁盘余量足够后让 r750-2 同时跑 C 和 D，但收益可能受 S3/网络带宽限制。

先确认三台机器看到同一个共享目录
===============================
在三台主机分别执行：
  R1=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_R1_v1_0_1_20261012
  sha256sum "$R1/inputs/F3V2_DELTA_CANDIDATES_A.tsv" "$R1/inputs/OLD_F3_ACCESSIONS.txt"
三台应完全一致。如果路径并非共享存储，先不要分布运行，改走 rsync/copy workflow。

安装
====
把本包解压到 ART_MathDiscovery_v1 根目录。然后：
  cd ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_DISTRIBUTED_ACQUISITION_R2_v1_0_2_20261013
  export F3V2_R1_REFERENCE=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_R1_v1_0_1_20261012
  ./tests/run_selftest.sh

运行示例
========
每台机器均 conda activate art_novelty，然后在同一个 R2 stage 目录：

as-4125:
  nohup ./02_acquire_stratum.sh A > logs/acquire_A.log 2>&1 &

r750-1:
  nohup ./02_acquire_stratum.sh B > logs/acquire_B.log 2>&1 &

r750-2:
  nohup ./02_acquire_stratum.sh D > logs/acquire_D.log 2>&1 &

D PASS 后 r750-2:
  nohup ./02_acquire_stratum.sh C > logs/acquire_C.log 2>&1 &

每层成功标志：
  PASS_F3V2_R2_STRATUM_A_COMPLETE
  PASS_F3V2_R2_STRATUM_B_COMPLETE
  PASS_F3V2_R2_STRATUM_C_COMPLETE
  PASS_F3V2_R2_STRATUM_D_COMPLETE

全部完成后仅在一台机器执行：
  ./02_merge_and_freeze.sh 2>&1 | tee logs/merge_freeze.log
应看到：
  PASS_F3V2_R2_DISTRIBUTED_54G_DOWNLOAD_AND_SOURCE_FREEZE_COMPLETE

之后才运行：
  ./03_run_54g_delta_scan.sh 12 4 2>&1 | tee logs/run_03_scan_console.log
  ./collect_stage05f3v2.sh

GPU
===
02 acquisition 不使用 GPU：S3 HEAD/下载 + zstd 解压 + FASTA 计数 + SHA256，主要是网络/CPU/I/O。
03 的 six-frame translation、HMMER hmmscan、MMseqs2 也都是 CPU 路线。H200 不会加速当前冻结实现；为了审计一致性不要临时换 GPU 算法。
