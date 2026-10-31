Stage05F3B4D-D shared-GPU resource hotfix, v1.0.0, 2026-10-11

用途：原 v1.0.0 的 D runner 将“GPU 必须无任何 compute process”作为保守资源门槛。该门槛不是科学协议的一部分。
当前 6 个结构扩展 target 仅 115-434 aa；原 idle-only D 尝试在任何 Boltz candidate 启动前即退出，因此可以在结果未知状态下冻结一个 resource-only hotfix。

新规则：
- target/sequence/Boltz2/MSA/seed/recycles/sampling/diffusion/mmCIF/R5C thresholds 全部不变；
- 仍然单 GPU、严格串行；
- 允许 GPU 上已有其他进程；
- auto 选择 memory.free 最大的 physical GPU；
- 启动时必须 >=60000 MiB free（可用环境变量 F3B4D_D_MIN_FREE_MB 调高，但不要调低以追结果）；
- OOM 仍标 technical/resource-censored，不当 negative。

运行：
  conda activate art_structure
  ./run_stage05f3b4d_d_boltz_shared.sh auto

Boltz 后继续使用原 Stage05F3B4D 目录的：
  conda activate art_novelty
  ./run_stage05f3b4d_d_foldseek.sh 32
  ./finalize_stage05f3b4d.sh
  ./collect_stage05f3b4d.sh
