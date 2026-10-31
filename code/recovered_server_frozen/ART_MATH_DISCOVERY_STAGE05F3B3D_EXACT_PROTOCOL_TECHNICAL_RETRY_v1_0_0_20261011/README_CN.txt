Stage05F3B3D — exact-protocol technical retry

目的：处理 F3B3C primary 32-family attempt 中的 10 个 Boltz TECHNICAL_FAILURE。

冻结规则：
- primary F3B3C 结果不覆盖、不重写。
- retry set = 全部且仅限 primary Boltz TECHNICAL_FAILURE，共 10 条。
- failure log diagnosis 仅描述，不参与 retry selection。
- 每条只做 1 次 exact-protocol technical retry。
- Boltz 参数完全不变：boltz2 / server-assisted MSA / seed 20261007 / recycling 3 / sampling 200 / diffusion 1 / --no_kernels / 1 GPU per job / mmCIF。
- 不截断序列，不使用 msa:empty，不调科学阈值。
- retry 成功后使用完全相同的 frozen R5C Foldseek 规则。

主结果与 technical rescue 必须在最终报告中分开记录。
