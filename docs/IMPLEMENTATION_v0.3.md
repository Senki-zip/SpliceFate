# SpliceFate v0.3 实现交付

2026-11-21 数据更新：GSE295932 已改为作者 `cluster` 整数 UMI 簇列，G 与 S 统一计数口径，重新选择 3,000 HVG；旧 count 加权版本保留但退出默认配置。细节、文库／细胞选择与散射过滤差异见 [数据转换开发记录](GSE295932_COUNT_SELECTION_20261121.md)。

2026-11-20 更新：官方 CytoTRACE 0.3.1 已复制至 `tools/CytoTRACE_official/`。默认对每个轨迹上下文使用完整基因计数评分，最高分细胞作为 DPT 起点；多阶段时在最早阶段内选择。所有事件共享同一份 S 校准时间。最新接口以 README 为准。

2026-11-20 后续更新：移除人工功能权重接口。新增 `composition.py`，以有公共锚点的 RBP 剪接状态和总表达共同拟合下游事件；从 WT 学习有符号权重，按冻结 WT 相对组成变化同步多轮传播。无固定时间 lag。蛋白结构域层仅作可选解释。

状态：核心代码实现完成，2026-11-16。依据用户最新指示，当前按可复用模型和独立合成数据验证，不继续根据当前数据集做模型取舍或启动完整训练。

后续澄清已实现：DPT 先排序，再仅对有可靠剪接信息的部分做局部小幅校准；共享/全局上下文也用 DPT，不切换 Palantir。没有剪接支持的细胞精确保留初始 DPT 值，不重新推断全局轨迹，不需要完整动力学。独立 `src/splicefate/intron_time.py` 提供计数方向与局部修正；完整 scVelo 适配仅保留为可选比较。最新验证为 17 项通过，新增测试确认部分支持时其他细胞原值不变；下文原 11 项为初次交付记录。

入口与完整输入/输出契约见项目 [README](../README.md)，配置见 `configs/mouse_v0.3.yaml`。

## 已实现

| 模块 | 入口 |
|---|---|
| 原始整数分子输入、HVG/归一化契约、身份/文件指纹 | `src/splicefate/io.py` |
| CausalBridge 配置 Leiden、小簇共享、统一 DPT 初始排序、断开状态 | `src/splicefate/trajectory.py` |
| 覆盖度加权 S 校准、全事件共享时间、可选 velocity 融合 | `src/splicefate/trajectory.py` |
| 内含子/外显子计数相位证据、仅支持细胞的局部 DPT 修正 | `src/splicefate/intron_time.py` |
| 原生 mm10 motif 扫描/CLIP、PPI 与复合体候选路径 | `src/splicefate/priors.py` |
| 联合条件主效应、样条、少量交互、DM 计数似然 | `src/splicefate/model.py` |
| 可判定事件–链映射、合法 simplex 重建、结构支持/冲突 | `src/splicefate/structure.py` |
| 单/联合剂量、冻结 WT 差分、同步多轮组成反馈与停止原因 | `src/splicefate/simulation.py` |
| prepare / fit / simulate / run / evaluate 命令，bootstrap 与逐事件审计 | `src/splicefate/cli.py` |
| 有真实 precursor/mature 输入时的 scVelo API 适配 | `src/splicefate/velocity.py` |

算法输入不要求固定细胞数或固定十个阶段；示例数据适配仍是小鼠原生 mm10。模块按输入类型启用，相位校准读取独立 E/I 计数。

## 验证

独立合成端到端输入包含原始 gene counts、事件分支、链结构、小型 DuckDB 先验和 FASTA。测试不读取项目真实数据，不访问网络。

11 项测试通过，涵盖：

- DM 目标函数解析梯度与有限差分一致。
- 非 log 单位上的剂量缩放、冻结其他因子、q=0 严格不变。
- 联合干预非加性、模型序列化一致、不涉及事件输入的干预严格无效应。
- 公共锚点缺失不当作排除分支、合法结构概率约束。
- 全部有支持的基因共同校准，所有事件共享一维时间，原 S 矩阵不变。
- 无 S 证据精确保留 G 轨迹、velocity QC 缺失阶段不改变时间区间定义。
- 没有真实 precursor/mature 量化来源时拒绝 velocity 输入。
- 组成项只使用实际观测，保留全部细胞，不填补缺失比例。
- WT 学习组成响应的同步多轮更新及 no-op 行为。
- 独立合成 prepare→fit→simulate 全流程。
- 外部评价必须匹配靶点、剂量、事件、状态、分支，不能混算。

以下为最初版本的历史验证记录（当时含基因块排除，现已移除）：`runs/implementation_validation/pytest.log`、`implementation_report.json`。scVelo 0.3.4 适配器另用独立 100 细胞、30 基因 U/S 数据跑通 recover_dynamics、velocity_graph、latent_time，按两个基因块排除并记录 16/14 个支持基因；见 `scvelo_adapter_check.json`。检查覆盖 API 和数值计算路径。

在用户要求“先不管当前数据集”之前，曾对现有输入完成一次 6 事件实现检查；该历史目录位于 `runs/splicefate_v0.3_implementation_check`。此后代码继续修订，没有继续用该数据集训练或按其结果选参；旧运行绑定当时源码，不能当作当前版本完整训练或复用缓存。

## 动力学模块和级联的明确边界

- scVelo 可选入口已实现并通过合成 API 检查；正式使用需独立验证的 spliced/unspliced 基因 UMI 与量化来源。
- Region Velocity 已检索并审查其长读长专用源码，尚未包装运行或重新量化 exon/intron regions。详见 [源码审查](LONGREAD_VELOCITY_SOURCE_REVIEW_20261116.md)。
- 直接 RBP 和伙伴模块通过其自身／伙伴 RNA 锚点生成事件候选。多轮传播由已拟合的 RBP 自身事件响应和下游组成系数连接，不要求人工蛋白功能注释。
- 可选 bootstrap 对完整 WT 事件模型集合重拟合并运行多轮组成传播，默认未开启时不输出区间。
- 没有实现 Δgene RNA 或命运概率；它们仍是 L1/L2 通过独立扰动验证后的后续模块。
