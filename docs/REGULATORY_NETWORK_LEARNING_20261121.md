# 全量单关系学习与质量筛选

更新日期：2026-11-25。

流程：先验候选 → 全部单关系学习与袋外评分 → 每事件最多 5 个输入 → 联合拟合 → 扰动模拟。学习顺序参考 [CausalBridge_en 的 fit_target / _wt_scores](https://github.com/Senki-zip/CausalBridge_en/blob/821c0314850fe12fde9d017d41c577bc4c7efd59/CausalBridge/saturating_hybrid.py)。

## 候选范围

`priors.targets: [ALL]` 从所有已表达剪接因子出发，按 `partner_hops: 1` 检索物理互作和同复合体伙伴。`extra_targets` 显式追加种子。检索跳数与 `simulation.max_rounds` 独立。

直接关系需要事件窗口内 motif 或 CLIP 锚点。伙伴关系需要伙伴在该窗口内的直接锚点，以及 recipient–partner 的物理互作或同复合体证据。来源证据分别保存于 `priors.duckdb` 的 `pair_evidence`、`event_anchors`、`event_candidates`；候选读取按事件进行。

## 完整与嵌套模型

单关系完整模型使用状态背景、一个 RBP 的总表达及有支持的自身局部剪接组成；预测目标基因时排除该基因的组成输入。嵌套参照保留相同状态列和源组成观测门控，移除表达、组成作用项。二者使用相同计数、细胞、固定浓度 `20.0` 和已冻结的 WT 设计基。

设计基只在目标事件有观测的 `good` 细胞上估计：表达尺度、样条节点、中心化和正交投影均固定。状态基去除完全共线列，时间在轨迹上下文内部标准化。唯一因子在观测细胞内无变异时记录 `constant_input_on_observed_cells`。分支数超过 8 的事件记录 `too_many_branches_unmodeled`。

## 袋外增益与稳定性

全部 WT 细胞用于最终模型拟合，无固定测试细胞集。网络筛选另外进行默认 10 次 `branch × stage × batch` 分层有放回重采样，各候选共用同一事件的样本。

每次完整、嵌套模型在袋内拟合，只在未抽中的袋外细胞评估：

- `oob_gain_b = mean(ll_full(C_oob, p_boot_oob) - ll_reduced(C_oob, p0_boot_oob))`，使用未加正则惩罚的 DM log likelihood。
- `agreement_b` 为全 WT 拟合的 `p_full-p_reduced` 与该次重拟合的响应向量在相同袋外行的正向余弦一致性。
- 增益超过 `min_oob_fit_gain` 时，该次稳定性取 `agreement_b`，否则取 0。
- 空袋外集和拟合失败记录原因，稳定性和正增益比例的分母仍为请求次数。

`oob_fit_gain` 是有效袋外拟合增益均值；`oob_positive_gain_fraction` 是超过增益门槛的次数除以请求次数；`stability` 是上述稳定性在请求次数上的均值。`fit_gain` 另保留全 WT 样本内增益作审计。尺度、样条、先验和伪时间在此步骤固定。

没有组成观测门控的候选共用同一事件的 state-only 全 WT 基线及重采样基线。带门控的候选按相同嵌套设计缓存基线；缓存不改变数值。

## 响应幅度与最终筛选

将 RBP 表达固定到观测 WT 均值、中心化组成对比置零，得到状态条件参照。完整概率减去该参照概率构成因子响应。各 `branch × stage × batch` 内按时间和 cell ID 排序，时间严格递增、两端均观测的相邻细胞累加响应差的总变差；分母使用全部时间递增相邻对，缺失位置不跨越。

入选资格：完整与嵌套拟合成功，`oob_fit_gain > 0`、`oob_positive_gain_fraction >= 0.6`、`stability >= 0.3`、响应幅度 `> 0`。阈值均在 `model.network_selection` 中配置。

合格候选的袋外增益和幅度各除以同事件合格候选最大值：

`selection_score = (0.45 × normalized_oob_gain + 0.45 × stability + 0.10 × normalized_response) / sum(weights)`

按综合分降序、factor ID 升序选择至多 `max_factors: 5` 个，允许 0～4 个。互作项只在入选因子内按两端质量得分排序，受 `max_interactions` 预算限制。本轮不做复合体小组救援。

## 并行、续跑与产物

`runtime.threads` 控制按事件多进程预算，数值 worker 内 BLAS/Numba 限为单线程。主进程独占网络写连接，事件关系与 `completed_events` 标记在同一事务提交。重启时先核对输入、代码、模型配置和事件范围指纹，跳过已完成事件，全部完成后统一排序。

- `regulatory_network.duckdb`：全部候选的拟合状态、紧凑系数、样本内与袋外增益、稳定性、重采样审计、排名和入选标记。
- `regulatory_network.tsv.gz`：仅结构化标识、数值和状态列，包含未入选及失败关系。
- `regulatory_network.summary.json`：候选、拟合、质量合格和入选数量。
- `models.duckdb` 与 `model_arrays.npz`：最终联合模型；预测和支持矩阵按事件偏移量读取。

`WT_weight` 是同一事件合格候选综合分的归一化值。扰动使用重新联合拟合的响应函数。`execution.bootstrap` 单独控制最终扰动区间，固定本次选择的输入网络，按事件收集成功的联合重拟合与级联结果。
