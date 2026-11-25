# SpliceFate v0.3 实现

更新日期：2026-11-25。改造依据：[第一轮完整运行方案](FULL_RUN_REDESIGN_PLAN_20261125.md)。

## 数据与主配置

主配置为 `configs/gse295932_stimulated_hg38.yaml`：GSM8961554 的 662 个刺激态 Nontarget WT 细胞独立转换、独立选择 Seurat v3 3,000 HVG，启用同条件 E/I 相位校准；Leiden 关闭。原合并数据和配置保留。

G/S 均采用作者 `cluster` 整数 UMI 簇数。刺激态转换得到 16,097 个基因、156,406 条观测剪接链、19,171 个事件和 40,756 个分支。完整 G 为 3,599,481 个簇；E/I 区域计数分别为 3,129,285 / 356,916，逐项复现 G 输入。

CytoTRACE 官方 0.3.1 对完整 G 评分确定 DPT 起点，再做有观测支持的 E/I 局部校准。全部 WT 细胞参与模型拟合，无固定测试细胞集，无 S 填补。

## 模块

| 模块 | 实现 |
|---|---|
| 按条件整数计数转换、Seurat v3 HVG | `scripts/convert_gse295932.py` |
| 同条件 E/I 区域计数与 G 逐项验证 | `scripts/quantify_gse295932_regions.py` |
| CytoTRACE、DPT 与 E/I 相位校准 | `cytotrace.py`、`trajectory.py`、`intron_time.py` |
| 独立 1 跳伙伴检索、motif/CLIP 锚点 | `priors.py` |
| OOB 质量评分、基线缓存、事件多进程与事务续跑 | `network.py` |
| 观测细胞设计基、状态列去共线、固定浓度 DM、联合支持域 | `model.py` |
| 观测组成特征、低计数平滑与同步级联 | `composition.py` |
| DuckDB、NPZ 内存映射与稀疏结构矩阵 | `storage.py` |
| 流式事件输出、三项贡献、状态与三种汇总 | `simulation.py` |
| 按基因并行的观测链重建 | `structure.py` |
| 同 WT 坐标的真实 KO 投影与独立评测 | `scripts/project_gse295932_ko_truth.py`、`cli.py` |
| 配置启动、R 包与官方源码预检 | `runner.py`、`cytotrace.py` |

源模块路径均相对 `src/splicefate/`。网络评分、默认门槛及续跑定义见[网络学习说明](REGULATORY_NETWORK_LEARNING_20261121.md)。

## 模型与输出约定

响应输入包括因子总表达、目标事件之外的自身观测组成、状态背景和少量受先验支持的交互。组成系数按下游事件从 WT 学习。`dm_concentration: 20.0`，事件最多 8 分支。

支持域使用本模型实际输入的表达、组成和上下文标准化时间。观测 WT 中的第 k 近邻距离分位数确定支持半径，k 默认 10。输出域外、低表达尾区和弱可辨识性标记。

组成传播以平滑参考 `(C + s × pooled)/(N+s)` 叠加相对 WT 的 log 概率变化，只更新有来源观测的细胞。零变化直接返回原输入。传播受轮数、收敛阈值和发散上限控制。

每个事件概率变化拆为表达剂量、直接组成缩放、后续级联三项，逐细胞加和等于总变化。细胞表默认只保存直接建模或实际级联改变的事件；全事件汇总持续输出，包含目标关系状态、网络质量排名和分子加权／观测细胞均值／全细胞均值。真实 KO 比较默认采用分子加权口径。

预处理保存 `priors_meta.json` 与 `priors.duckdb`。最终模型保存 `models.duckdb` 与连续数组 `model_arrays.npz`；剪接链保存 `structures_meta.json` 与稀疏 `structure_counts.npz`。各阶段以 SHA256 绑定输入、源码和产物。

## 外部 KO

投影脚本仅使用 WT 已发现的公共锚点和分支。官方 barcode 注释区分 11 个剪接因子 KO、BUB3 外部覆盖审计及 VIM 阴性对照；真实 KO 不参与网络学习或选参。WT 重投影结果须与训练分支计数逐项相同。外部真值文件记录分支分子数、WT/KO 观测细胞数与过滤状态，实际 KO 对应 `strength: 1.0`。

## 验证

本地测试覆盖固定设计的袋外行隔离、噪声候选增益、缓存恒等、串并行一致性、事务续跑与指纹拒绝、DuckDB/NPZ/CSR 读写、观测细胞设计基、零剂量恒等、三项贡献加和、支持域、状态与汇总口径，以及 R 依赖检查。测试文件仅保留在本地。

真实数据验证采用刺激态全部 662 个 WT 细胞与少量事件，执行 prepare→fit→11 KO×3 剂量 simulate→evaluate。完整事件运行由主配置单独启动。

本轮验证完成：111 项测试、662 WT×20 事件、11 KO×3 剂量及真实中断续跑通过；完整统计见[验证记录](FULL_RUN_REDESIGN_VALIDATION_20261125.md)。
