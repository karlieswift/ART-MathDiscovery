# Stage05D4-R5 v1.0.1 implementation hotfix: 4 independent single-GPU shards

Reason: the original `--devices 4` Boltz/PyTorch-Lightning run entered multi-process NCCL/DDP. One rank exited and the remaining ranks emitted TCPStore/NCCL broken-pipe warnings. The one-GPU smoke test had already passed.

This hotfix avoids distributed NCCL entirely. It partitions the frozen 43 P1 YAML inputs deterministically (sorted round-robin) into four shards and runs four independent Boltz-2 processes, one per H200 GPU.

Scientific settings are unchanged: Boltz-2, MSA server, seed 20261007, recycling=3, sampling_steps=200, diffusion_samples=1, mmCIF output, `--no_kernels`. No DNA score/rank, candidate universe, P1/P2 rule, or biological threshold is changed.

The failed multi-GPU output is preserved by renaming `boltz_out` to `boltz_out_ddp_failed_<timestamp>` before canonical output assembly.
