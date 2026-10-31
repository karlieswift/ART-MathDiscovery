J2：外部评价与匹配对照的执行包
版本 J2-1.0.0，2026-10-22

一、这次实际交付的边界

本目录提供可执行的冻结、原算法评分、评分锁定和后验评价程序；不是把已有13例、UG27、919、6.12G或54G重新命名成外部盲测。正式ART外部队列尚需独立作者真值和预先匹配的对照。空队列会返回 HOLD_EXTERNAL_ADMISSION（退出码3），不会输出虚假的科学PASS。

另支持 secondary_structured_rna 模式。其他结构RNA的作者注释可检验数学结构定位的迁移；它们既不是ART阳性，也不是ART真阴性。这些一般结构RNA位点未必含RT，所以RT路线在这些位点无命中不能解释成ART检测性能差；这一次要队列不能支持“Arm M的ART检测能力优于RT方法”的结论。正式ART方法比较仍需独立的primary_art主队列。旧项目曾提及Minerva，未能证明过去未见的新来源保留 EXPOSURE_UNCERTAIN / SECONDARY_ONLY；不能声称新的盲外部验证。

已知曝光索引见 reference/EXCLUSION_LEDGER.tsv 和来源hash表。它包括已提供的正式/开发/回顾性基准、UG27、919、历史和54G源表的accession、原生contig标识、可用DNA/反向互补hash、蛋白hash及操作性family ID。索引中的“未找到”不构成从未见过的证据；需要独立的来源/相似性与家族分组审核。蛋白hash、局部DNAhash、完整基因组hash分别保留，不能混为一个完整基因组的证据。

二、算法没有换成简化版

vendor/ 保存按字节复制的原始 arm_m、clean_two_arm_replay 和 formal_baseline 代码。reference/VENDORED_SOURCE_SHA256.tsv记录每个文件；freeze前全部校验。默认 exact_two_arm 模式用新协调器直接调用原 make_windows、init_acnv/acnv_one、score_block、empirical_upper_q、recurrence_scores，再调用原 build_two_arm_board 构建完整 ACNV+generic recurrence 两分支 OR 区域及原排序；原 main() 内联的 family 校准运算保持相同操作和顺序。旧 IRR 分支完全不参与该正式方法，默认不执行。可选 --execution-mode full_legacy 仍执行原完整三分支扫描，随后重建正式两分支 board；不使用依赖 IRR 的旧三家族 union。详见 docs/EXACT_TWO_ARM_DEPENDENCY_AUDIT.txt。没有调阈值、使用truth裁剪窗口或筛选top-K幸存者。原引擎内部输出保留其历史Stage03A命名和旧状态文字；它们不能重新认证新数据的盲态，科学身份以外层analysis_kind、admission和EVALUATION_STATUS为准。

三个主要方法是 ARM_M_FULL_FROZEN、PRODIGAL_RVT_PLUS_FROZEN_REC、SIXFRAME_RVT_PLUS_FROZEN_REC。RVT-only两个路线作为次要比较。ACNV-only/recurrence-only是既有研究已完成的消融，不在本包冒称新补的缺失工作。RT对照调用原score_case、RT-locus合并和原rt_orders。Prodigal使用原meta设置；不强加新的translation table。sixframe沿用原ORF最小45aa规则。RVT模型full E和domain i-E阈值均1e-4，recurrence q为0.005，局部邻域±10kb。

三、输入与资源

每行一个完整原生contig，所有方法共用同一FASTA。版本1只接受单条、ACGT-only、至少4096bp的记录，不拼接多contig，不静默删N、不把已知truth附近窗口当完整基因组。多contig来源要在评分前列出完整抽样/纳入规则；按native parent、source study、family和匹配组的连通依赖单位分组，不能只取有truth的contig后称全基因组评价。

实际评分需要Python、numpy、pandas、tqdm、Prodigal、HMMER的hmmscan和hmmpress；不用GPU或结构预测。默认精确两分支模式的临时float32 memmap字节数为 oriented_windows × 5468 ×4。已提供6,722,539 bp全基因组有827,154个oriented windows，ACNV memmap为18,091,512,288 B（约18.1 GB），另需序列/特征数组、评分工作区、结果和进程内存。可选full_legacy模式还计算旧IRR，memmap为29,903,271,408 B（约29.9 GB）。规划RAM启发值约28.3 GB，不能拿它冒充实测峰值；请预留额外磁盘和足够RAM。admit的inventory给出每个输入的确切窗口数和memmap大小；不将它当峰值内存测量。CPU/超时在freeze时固定；失败样本不得在看结果后无记录删除。

必须使用原RVT-All.hmm，SHA256：
74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e
默认定位原项目Stage05D1、R5B、F3B中的已存在文件。若缺失，fetch-rvt可从官方myRT源恢复同一字节文件，下载后必须精确hash匹配；不是改用当前Pfam或新模型。模型沿用上游许可，本包不替模型重新授权。

四、命令（输出均自行指定到本项目内）

ART_PROJECT=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
J2_KIT=/实际解压路径/j2
J2_OUT="$ART_PROJECT/outputs/J2_20261022"
mkdir -p "$J2_OUT"

1. 环境检查
bash "$J2_KIT/run_j2.sh" preflight --project-root "$ART_PROJECT" --output "$J2_OUT/preflight"

如Prodigal在其他已有conda环境，加入 --prodigal /完整路径/prodigal。同理支持 --hmmscan、--hmmpress 和 --rvt；无需改变原环境PATH。

2. 只在原HMM缺失时精确恢复
bash "$J2_KIT/run_j2.sh" fetch-rvt --output "$J2_OUT/reference/RVT-All.hmm"

3. 检查新的manifest和预先固定的独立真值
bash "$J2_KIT/run_j2.sh" admit --candidates /路径/candidates.tsv --truth /路径/truth.tsv --output "$J2_OUT/admission"

输入字段见 config/CANDIDATES.template.tsv、TRUTH.template.tsv 和 MANIFEST_SCHEMA.json。所有作者证据和leakage review均带本地文件hash。真值必须来自独立作者注释/认证，不能根据本研究detector的结果定义。primary_art的每个阳性match_group必须含独立定义的hard_neighbor或certified_negative；随机背景不能自动变成阴性。background_unknown只报告候选负担。

secondary结构RNA使用同样的admit命令并添加 --analysis-kind secondary_structured_rna。标签为 other_structured_rna、author_defined_rna_interval、INDEPENDENT_AUTHOR_ANNOTATION；没有RT真值不会阻挡这个次要迁移分析。同一genome的多个作者位点分别target_id，分母仍显式报告native genome数和独立依赖单位数。

4. 冻结（run目录必须尚不存在）
bash "$J2_KIT/run_j2.sh" freeze --project-root "$ART_PROJECT" --candidates /路径/candidates.tsv --truth /路径/truth.tsv --run "$J2_OUT/run_001" --cpu 16

如是secondary，freeze也必须添加 --analysis-kind secondary_structured_rna。如需同一环境外的工具或恢复的HMM，添加对应--prodigal/--hmmscan/--hmmpress/--rvt参数。freeze固定execution mode、协调器SHA，复制完整输入、原始源码和所有参数，重建并hash锁定HMM索引，记录工具binary hash、版本和Python依赖版本。truth、原始身份键、类别和来源审核放在sealed/；scoring manifest只含opaque ID和完整输入信息。

5. 真实运行全部方法并锁定输出
bash "$J2_KIT/run_j2.sh" score --run "$J2_OUT/run_001"

评分函数不会打开sealed/，原始输出和每个case-method成功/失败状态全部保留，随后锁定hash。两条RT路线独立执行：Arm M失败仍运行两条RT-only；只有真实recurrence输入缺失的RT+rec记dependency failure，不能把数学引擎失败算成所有对照失败。每阶段保留wallclock、child CPU时间和Linux wait4或GNU time -v测量的最大child RSS；RSS不是多worker同时占用内存的总和。中途退出尚未写出完整SCORE_FREEZE.json时，不能开truth。默认不设置任意运行时截断；如需要上限，必须在freeze时用--timeout-seconds预先指定（每个引擎阶段上限），并如实保留超时失败。

6. 全部方法的评分尝试锁定以后才能评价
bash "$J2_KIT/run_j2.sh" evaluate --run "$J2_OUT/run_001"

evaluate先验证每个case-method的完整状态和所有输出hash，再验证并读取sealed真值。缺方法、改结果、改参数或改truth都会HOLD，不自动补成成功。失败方法仍在完整阳性分母中按未恢复计；失败样本的候选负担为缺失，不错误写成零。软件的文件隔离是计算流程防泄漏，不代表作者在现实中不知道已有注释。

五、报告的端点与解释

TARGET_METRICS/TARGET_SUMMARY报告预先定义目标区间的Hit@1/3/5、recall、MRR。ART主目标是作者独立定义的system interval，阳性输入还必须同时提供该system内独立的array和RT坐标，组件分别在COMPONENT_LOCALIZATION报告；不能把system-span overlap写成精确ncRNA或功能验证。每个额外target执行与第一行同等的scope、坐标和证据检查。其他结构RNA模式只评价作者RNA位点，不混入ART阳性统计。已暴露family会明确记录FAMILY_EXPOSED；本包不会仅凭accession不同认证family-disjoint generalization。

ALL_SAMPLE_CANDIDATE_BURDEN保留所有阳性、hard-neighbor、certified-negative和unknown背景的候选数/Mb、提名区间占基因组比例、RT±10kb暴露比例。这里不会把无真值的随机背景解释成false-positive/precision。PER_GENOME_TARGET_METRICS和GENOME_AGGREGATE_SUMMARY保留基因组分母，避免把同一基因组几十个位点当几十个独立外部样本。

bootstrap按freeze前固定的leakage_unit做同一重采样下的方法配对比较，seed=20261022、2000次。独立阳性依赖单位少于5只给描述结果，不生成虚假的稳定置信区间；达到最低值也不证明样本具有代表性。所有方法同一原生输入，但候选暴露宽度可能不同，必须同时看区域/邻域暴露比例。

六、真实历史集成验证与软件防护测试

bash "$J2_KIT/run_j2.sh" historical --output "$J2_OUT/historical_SA1" --cpu 16 --rvt /路径/RVT-All.hmm --prodigal /路径/prodigal --hmmscan /路径/hmmscan --hmmpress /路径/hmmpress

该命令使用包内完整SA1原生基因组260727bp，默认真正重算完整原始ACNV和recurrence；可选full_legacy另外执行不参与正式方法的旧IRR，构建原两分支排名，并运行实际Prodigal/HMMER和sixframe，随后逐项核对全部32,050个oriented窗口键/顺序、16,025个genomic窗口和全部非IRR数值、78个canonical区域全部字段、160条canonical proposal及原RT位点/数值。它明确是historical regression，不是新的外部性能。编译得到同版本工具不保证与原binary逐字节一致，因此报告binary hash并比较实际数值，差异不会假装PASS。

bash "$J2_KIT/run_j2.sh" test --output "$J2_OUT/SOFTWARE_GUARD_TESTS.json"

此命令只检查缺方法、hash变化、多contig、N、失败分母、小样本伪重复、空外部表等软件边界；临时合成fixtures不是科学分析结果。真实历史CPU回归由上一个命令完成。

七、本包附带的实际验证结果

validation/deliverable/HISTORICAL_REGRESSION_STATUS.json 记录真实完整SA1回归：32,050个oriented窗口、16,025个genomic窗口、35个非IRR列共560,875单元格全部exact一致；78个canonical区域1,638字段和160条proposal全部一致。实际Prodigal和sixframe/HMMER给出2个旧RT位点，12个证据数值一致，6项独立/缺rec路线检查一致。原始扫描wall为164.374秒（8 workers），最大child RSS为1,462,608 KiB；这不是多worker总峰值或6.72Mb外部样本的测量。

同目录保留所有3次先前失败/停止尝试的原因；新成功的完整精确两分支运行不抹去这些记录。13个软件防护测试通过，但它们不是新外部科学结果。已提供的secondary数据仍需服务器完整评分；正式primary_art队列当前无新合格样本。
