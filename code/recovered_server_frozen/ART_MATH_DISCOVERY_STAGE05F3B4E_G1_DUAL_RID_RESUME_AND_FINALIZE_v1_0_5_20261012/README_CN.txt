Stage05F3B4E G1 R5：双 RID 恢复 + 自动 finalization。
用途：停止旧 R4 本地 WAITING 轮询后，只复用已经提交的三个 NCBI RID：
RT5 CR73H1ST014；partner5 CSE25BRU016；architecture5 CSRWCA1M016。
不会提交任何新 RID，不改任何科学参数。partner 与 architecture 并行等待/取回；完成后自动 summarize_global，复制 G1 canonical outputs 到 base F3B4E，随后运行已有 finalize_stage05f3b4e.sh 与 collect_stage05f3b4e.sh。
若 NCBI 仍 WAITING，保持运行即可；若本地中断，重新运行同一命令继续使用相同 RID。
