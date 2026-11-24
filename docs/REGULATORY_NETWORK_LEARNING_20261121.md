# 全量单关系学习与每事件前 5 个输入

更新日期：2026-11-24。

流程为：先验候选网络 → 全部单关系拟合与 WT 重采样 → 增益／稳定性／幅度综合评分 → 每事件最多 5 个输入 → 联合拟合 → 扰动模拟。

## 学习顺序

参考 [CausalBridge_en](https://github.com/Senki-zip/CausalBridge_en/blob/821c0314850fe12fde9d017d41c577bc4c7efd59/CausalBridge/saturating_hybrid.py) 的 `fit_target`、`_wt_scores`：先拟合所有单输入关系，再选择组合输入。SpliceFate 在响应幅度之外加入条件计数拟合增益与 WT 重采样稳定性，默认分别占 45%、45%、10%。这组比例为可配置的初始工程权重。

SpliceFate 的第一遍处理遍历全部候选 RBP–事件关系。每条关系使用全部 WT 细胞，计数目标由该事件有观测的细胞提供；输入保留状态背景、单个 RBP 总表达及有支持的自身局部组成。自身组成输入使用原有公共锚点规则，目标基因的同一组成信息从自身预测输入中排除。

所有单关系保存完成后进入第二遍处理，按学习得分为每个事件选择最多 5 个因子，联合拟合原有条件 DM 模型。先验的 seed、CLIP/motif/伙伴类型及 ID 顺序保存在证据中，入选顺序由 WT 学习得分决定；完全相同的得分按稳定因子 ID 排序。

`priors.targets: [ALL]` 定义全表达候选网络的起始集合。`simulation.target_sets` 定义后续查询靶点，两者分别用于网络构建和扰动。

## 拟合增益与稳定性

每条关系先拟合原有条件 DM 响应，再拟合嵌套参照：保留同一状态设计列、相同的源组成观测门控列，移除 RBP 表达和组成效应。两者使用相同的事件观测细胞、计数和固定 DM concentration。拟合增益为每个观测细胞的未加正则惩罚 DM log likelihood 差的均值，省略项均为相同的参数无关计数常数。

稳定性使用 WT 内 `branch × stage × batch` 分层有放回重采样，默认 10 次。每次保留原层的观测细胞数，在固定 WT 特征尺度、样条基和组成门控下，重新拟合完整关系和嵌套参照的全部系数。将两模型重新求值到共同的原 WT 观测细胞，得到 bootstrap 拟合增益和 `p_full-p_reference` 响应向量。每次分数为：增益超过 `min_fit_gain` 时取该向量与全 WT 响应向量的正向余弦一致性，否则取 0；数值拟合失败也贡献 0。稳定性取全部请求重采样的平均，范围为 0～1。

该指标衡量固定 WT 表示下的拟合重现性，所有 WT 细胞用于最终拟合与比较。网络拓扑、先验和伪时间在此步骤固定。分层重采样在候选间共用，随机种子由配置和稳定事件 ID 决定。

## WT 响应幅度

对关系 r→e，先计算每个细胞的单关系分支概率 `p_i`。然后固定 r 的 log1p 归一化表达输入为该事件观测细胞的 WT 均值，将 r 的中心化组成对比设为零，保留细胞原有状态和观测门控，得到参照概率 `p_ref_i`。差值 `d_i = p_i - p_ref_i` 表示拟合的因子响应。

在每个 `branch × stage × batch` 中，按共享伪时间、cell ID 排序。伪时间递增的相邻细胞形成比较对；只在两端均有事件计数时累加 `0.5 × sum_k |d_jk-d_ik|`，最终除以这些上下文中全部可比较相邻对的数量。缺失事件观测保留原位置，分母在同一事件的所有候选间一致。

此处将 CausalBridge 标量响应的相邻绝对变化扩展为多分支概率的总变差，并通过配对参照分离状态背景。相同时间、跨上下文、跨阶段和跨批次的细胞分别处理。

## 综合评分与筛选

拟合成功、嵌套参照可用、拟合增益超过 `min_fit_gain`、响应幅度超过 `min_response_score` 且稳定性大于 0 的关系参与排序；两个阈值默认均为 0。在同一事件合格候选内，将增益和响应幅度分别除以对应最大值。稳定性直接使用其 0～1 值。

`selection_score = 0.45 × normalized_fit_gain + 0.45 × stability + 0.10 × normalized_response`

权重写入 `model.network_selection.weights`；按配置权重总和归一化。联合模型保留综合得分前 `max_factors` 个，默认最多 5 个。已有互作项在入选因子内按两端综合得分排序，预算由 `max_interactions` 控制。未满足质量条件的关系保留在完整网络记录中，三项原始评分分别输出。

配置方法为 `all_single_relations_then_weighted_quality`；`bootstrap_repeats` 和 `bootstrap_seed` 控制网络选择重采样，独立于 `execution.bootstrap` 的最终扰动区间重采样。修改评分配置时使用新输出目录重新拟合。

## 输出

- `regulatory_network.duckdb`：`single_relations` 保存每条候选的事件、因子、拟合状态、响应幅度、支持比较对数、单关系曲线、来源证据、入选排名和标记；质量字段包括 `quality_status`、`fit_gain`、`stability`、`bootstrap_successes`、`bootstrap_positive_gain_fraction`、`bootstrap_effect_agreement`、`bootstrap_gain_mean`、`normalized_fit_gain`、`normalized_response`、`selection_score`。`quality_json` 另保存基线 likelihood，`network_metadata` 保存评分权重及重采样配置。
- `regulatory_network.tsv.gz`：导出全部关系及筛选结果，含未入选关系与拟合失败状态。
- `regulatory_network.summary.json`：候选、成功拟合、正得分、入选关系和状态计数。
- `models.json`：最多 5 因子的最终联合响应模型，用于扰动和组成级联。

单关系曲线保存系数、设计基和拟合诊断；细胞级输入统一从准备阶段产物读取。`WT_weight` 是同事件合格关系综合得分的归一化值，用于网络解释，最终扰动使用重新联合拟合的响应函数。`execution.bootstrap` 固定本次 WT 选择的网络输入，并重拟合联合模型和组成级联。

全量执行配置为 `execution.event_limit: null`。小规模检查可限制事件数量；每个被检查事件的全部候选关系仍先完成学习。线程、归一化、HVG、模型正则化及 E/I 输入设置沿用原配置。
