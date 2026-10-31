Stage05F3B4E Public ART8 R1
原因：base runner 在任何 ART8 BLAST 结果生成前因 missing unzip 退出。
本 hotfix 只把 unzip 命令替换为 Python 标准库 zipfile；不改变 5 targets、8 public ART genomes、datasets request、BLAST 参数或 downstream thresholds。

Run:
  conda activate art_novelty
  ./tests/run_selftest.sh
  ./run_stage05f3b4e_public_art8_r1.sh 2>&1 | tee logs/run_public_art8_r1_console.log
