F3B4E G1-R1 result-blind remote transport hotfix.

用途：原 G1 batched remote BLASTP 在任何结果产生前两次 transport failure。
本包不改变 query、database scope、evalue、max_target_seqs、SEG/DUST 或 novelty thresholds。
只改为每条 query 单独提交，并用 6 次退避重试；若本机恰有同名 nr/nt BLAST DB，则允许相同参数 local fallback。
任何 transport failure 都保持 technical-censored，绝不当 no-hit。

依赖：原 F3B4E prepare 已成功。
