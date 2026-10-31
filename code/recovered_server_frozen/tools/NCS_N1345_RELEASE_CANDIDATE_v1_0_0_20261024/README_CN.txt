ART MathDiscovery NCS N1/N3/N4/N5 v1.0.0 - 2026-10-24

这是科研代码与源数据候选包（不是论文或投稿信）。

实际核验：
N1=HOLD 目前没有可认定为新的、家族独立、有确切ART阵列+RT坐标的基准真值；不制造结果。
N3=PASS 6基因组、3q阈值、18组后验区域/排名敏感性分析。
N4=PASS 67个结构x32参考，2144条结构对齐的归一化敏感性与质量审计；2个原历史失败，人工审核未完成。
N5=PASS_LOCAL 8个数据输入哈希、7项自动测试、7项输出表逐字节复算；没有实际创建Github或DOI。

解压后运行：
  conda activate art_novelty
  python -m unittest discover -s tests -v
  bash run_all.sh [输出文件夹]
  python src/visualize_and_export.py --root . --results [输出文件夹] --out [图/表输出文件夹]

注意：研究源表中的源代码路径已经净化；原始J05工作簿仍在作者的历史内部审计文件中。
本包无生化实验、无新增独立ART检测成绩，不改变现有论文、图稿及冻结主阈值。

全窗口矩阵级重放（需 pandas、numpy）：
  bash run_full_n3_matrix_replay.sh [输出目录]
  正常输出 PASS_FULL_SIX_MATRIX_RECOMBINATION_TWO_OUTPUTS_BITWISE_MATCH

