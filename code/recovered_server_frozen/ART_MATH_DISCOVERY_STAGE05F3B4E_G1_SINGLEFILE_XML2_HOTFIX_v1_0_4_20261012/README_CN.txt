Stage05F3B4E G1 R4：READY RID single-file XML2 retrieval hotfix。

R3 已证明 RID CR73H1ST014 为 READY，但反复 READY。原因是 R3 请求 FORMAT_TYPE=XML2（NCBI UI/URL API 中为 multiple-file XML2）后，却按 single-file XML tags 判断返回值，因此把“结果已完成但取回格式不匹配”误当作格式暂未就绪。

R4 仅将 Get 的结果序列化从 XML2 改为 XML2_S（single-file XML2），不重新提交已有 READY RID，不改变任何 BLAST search 参数、query、database 或 novelty threshold。
