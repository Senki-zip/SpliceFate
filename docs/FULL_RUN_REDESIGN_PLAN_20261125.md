# SpliceFate 第一轮完整运行改动设计方案

日期：2026-11-25

适用基线：提交 `b9e1f03`（v0.3）

目标：在保留核心科学约束（全部 WT 细胞参与拟合、不划分测试集、无 S 填补、整数 UMI 簇口径、CytoTRACE 0.3.1 定向 DPT、先全量单关系学习再筛选前 5 联合拟合）的前提下，解决全量运行的计算、内存/存储、统计评分偏置、建模细节及外部 KO 评测对齐问题，为第一轮完整运行做好准备。

---

## 1. 定稿决策总览

| 编号 | 模块 | 最终决策 |
|---|---|---|
| 1 | 候选范围与计算扩展 | `priors.partner_hops: 1` 与 `simulation.max_rounds` 解耦；缓存事件级 state-only 基线与重采样；24 核按事件多进程并行；`regulatory_network.duckdb` 按事件增量落盘与断点续跑 |
| 2 | 存储与模拟输出防 OOM | `priors`/`models`/`structures` 拆离大 JSON，改为 DuckDB + 二进制 `.npz`；`regulatory_network.tsv.gz` 仅导结构化数值列；模拟表流式分块写盘，默认只写受影响事件的细胞级行（`write_unchanged_cells: false`），`event_stage_summary.tsv` 始终全量输出；`structure_predictions` 按基因并行 |
| 3 | 质量评分与筛选门槛 | 在已有 10 次 WT 分层有放回重采样上启用袋外（OOB）增益与 OOB 稳定性，沿用 45% / 45% / 10% 权重；门槛 `min_oob_fit_gain: 0.0`、`min_positive_gain_fraction: 0.6`、`min_stability: 0.3`，未达标事件允许入选 `<5` 个因子；第一轮不做复合体小组救援 |
| 4 | 数据集拆分与核心建模 | 拆开 `stimulated` 与 `unstimulated` 独立建模，第一轮主数据集采用 `stimulated`（GSM8961554，662 WT 细胞）；通用代码清理 `state_basis` 完全共线列并按上下文内部标准化时间；设计基改在观测细胞 `good` 上构造；组成级联低计数平滑；三项贡献拆分 + 联合支持域修复 + 悬空配置清理；固定 `dm_concentration: 20.0` 且 `max_event_branches: 8` |
| 5 | 输出状态、评测与运行配置 | `event_stage_summary.tsv` 增加 `target_relation_status` 与三种聚合口径；新增真实 KO 同坐标投影脚本对接 `evaluate`；主配置开启 E/I 相位校准（`stimulated` 专属 `regions/`）；默认批量模拟 `stimulated` 中 11 个剪接因子 KO（`BUB3`/`VIM` 单列外部审计）；修复全量 `execution.bootstrap` 与 `--check` R 环境校验 |

---

## 2. 分模块详细设计

### 2.1 数据集拆分与 E/I 区域量化（`scripts/convert_gse295932.py`、`scripts/quantify_gse295932_regions.py`）

1. **按条件独立转换**
   - `scripts/convert_gse295932.py` 新增 `--condition {stimulated,unstimulated,both}` 参数：
     - `stimulated`：仅导入 `GSM8961554_FT_stim_sc_iso_count.Nontarget.txt.gz`（662 个 WT 细胞），独立计算完整文库、Seurat v3 3,000 HVG、剪接链与局部公共锚点事件。
     - `unstimulated`：仅导入 `GSM8961552_FT_sc_iso_count.Nontarget.txt.gz`（387 个 WT 细胞），按需可独立转换与运行。
   - 默认输出目录按条件隔离，例如：
     - `data/processed/GSE295932_Nontarget_stimulated_hg38_cluster/`
     - `data/processed/GSE295932_Nontarget_unstimulated_hg38_cluster/`
   - 保留原合并目录不被覆写。

2. **单条件专属 E/I 区域量化**
   - `scripts/quantify_gse295932_regions.py` 同步支持 `--condition` 与 `--processed-dir`，对 `stimulated` 的 662 个 WT 细胞生成配套 `regions/`（`exon_counts.npz`、`intron_counts.npz`、`boundary_counts.npz`、`retained_counts.npz`、`intron_without_full_retention_counts.npz`、`cells.tsv`、`genes.tsv`、`quantification.json`），校验总 G 计数与 `genes/gene_counts.npz` 逐项一致。
   - 主配置中 `input.intron_counts` 指向该 `regions/` 目录，`trajectory.S_calibration: intron_phase` 正式生效。

---

### 2.2 候选先验范围与存储重构（`src/splicefate/priors.py`、`src/splicefate/cli.py`）

1. **解耦 PPI 扩展跳数**
   - 新增配置项 `priors.partner_hops`（非负整数，默认 `1`），不再借用 `simulation.max_rounds`。
   - 新增可选配置项 `priors.extra_targets: []`（用于显式追加额外查询靶点到候选宇宙，不改变默认 1 跳范围）。
   - `targets: [ALL]` + `partner_hops: 1` 下：从全部已表达剪接因子种子（0 跳）向外仅扩展 1 跳直接物理互作（STRING physical ≥ `min_experiment_score` / IntAct）或同复合体（CORUM / ComplexPortal）伙伴，剔除 2–3 跳远端非剪接蛋白。
   - 间接伙伴路径仍保持 1 跳传递语义：仅当 `partner` 在该事件窗口内有直接 motif/CLIP 锚点时，`recipient` 才通过 `(recipient, partner)` 边成为 `indirect_partner_candidate`。

2. **先验存储去单文件大 JSON**
   - 取消将全量重复证据写入单个 `priors.json`，拆分为：
     - `priors_meta.json`：轻量元数据（`factor_ids`、`gene_ids`、`symbols`、`seed_factor_ids`、`factor_retrieval_depth`、`eligible_events`、`isoform_states`、`composition_features`、`isoform_state_policy`、扫描与配置摘要）。
     - `priors.duckdb`：规范化关系表，避免重复拷贝 `relation_evidence`：
       - `pair_evidence(pair_id, factor_a, factor_b, max_rank, evidence_json)`
       - `event_anchors(event_index, factor_id, anchor_count, has_clip, anchors_json)`
       - `event_candidates(event_index, factor_id, gene_id, symbol, route, priority, anchor_count, partner_ids_json, path_count)`（按 `event_index` 建索引）
   - `prepare` 产物指纹 `provenance.json` 记录 `priors_meta.json` 和 `priors.duckdb` 的 SHA256。
   - 保留内存读取适配函数，供单测和小规模脚本按 `event_index` 批量查询，不再整库加载进 RAM。

---

### 2.3 网络学习、OOB 质量评分、并行与断点续跑（`src/splicefate/network.py`）

1. **袋外（OOB）质量评分与自由度偏置消除**
   - 保持“全部 WT 细胞参与最终拟合、不划分测试集”原则不变。
   - 在每个事件的 `B = bootstrap_repeats`（默认 `10`）次 `branch × stage × batch` 分层有放回重采样中：
     - 每次抽样天然产生约 63.2% 的袋内索引 `in_bag` 与约 36.8% 的袋外索引 `oob`（排除 `oob` 为空的极端抽样）。
     - `full` 与嵌套 `reduced` 模型在 `in_bag` 上拟合，在 **`oob` 细胞上评估**：
       - `oob_gain_b = mean(ll_full(C_oob, p_boot_oob, kappa) - ll_null(C_oob, p0_boot_oob, kappa))`
       - `oob_agreement_b = max(0, cos((p_wt - p0_wt)_oob, (p_boot - p0_boot)_oob))`（仅当 `oob_gain_b > min_oob_fit_gain` 时计入，否则贡献 `0`）
     - 汇总指标：
       - `fit_gain`：保留全 WT 样本内未惩罚 DM 对数似然差用于审计对照。
       - `oob_fit_gain`：所有有效 OOB 重采样中 `oob_gain_b` 的均值（若拟合失败按失败记录，不伪造正增益）。
       - `oob_positive_gain_fraction`：`oob_gain_b > min_oob_fit_gain` 的比例（分母为请求重采样次数 `B`）。
       - `stability`：`oob_agreement_b` 在请求重采样次数 `B` 上的均值（范围 `[0, 1]`）。
   - **入选资格与综合评分**：
     - 资格条件（全部写入 `model.network_selection` 配置）：
       - `status == 'ok'` 且 `quality_status == 'ok'`
       - `oob_fit_gain > min_oob_fit_gain`（默认 `0.0`）
       - `oob_positive_gain_fraction >= min_positive_gain_fraction`（默认 `0.6`）
       - `stability >= min_stability`（默认 `0.3`）
       - `response_score > min_response_score`（默认 `0.0`）
     - 事件内合格候选归一化：
       - `normalized_fit_gain = oob_fit_gain / max(oob_fit_gain)`
       - `normalized_response = response_score / max(response_score)`
       - `selection_score = (0.45 * normalized_fit_gain + 0.45 * stability + 0.10 * normalized_response) / sum(weights)`
     - 按 `selection_score DESC, factor_id ASC` 保留最多 `max_factors`（默认 `5`）个因子；**未达标事件允许入选 `0~4` 个因子**，避免噪声事件硬塞满 5 个。
     - 第一轮不加入复合体小组联合增益救援。

2. **事件级 state-only 基线缓存**
   - 同一事件下，若候选因子无 `composition_masks`，其全 WT 嵌套基线 `baseline` 以及 10 次重采样中的 `reduced` 基线完全相同。
   - 改为在每个事件循环开头预先拟合一次 state-only 基线及其 `B` 次重采样结果并缓存；仅对携带 `composition_masks` 的少量候选单独拟合含门控列的嵌套基线。
   - 数学结果严格等价，单关系阶段 L-BFGS-B 计算量约减少一半。

3. **按事件多进程并行与断点续跑**
   - 使用 `ProcessPoolExecutor(max_workers=runtime.threads)` 按事件分批并行执行单关系拟合与 OOB 评分；worker 进程启动时强制设置 `OMP_NUM_THREADS=1`、`OPENBLAS_NUM_THREADS=1`、`MKL_NUM_THREADS=1`、`NUMBA_NUM_THREADS=1`，避免 BLAS 线程争抢。
   - 主进程独占 DuckDB 写连接，按完成事件批次事务写入 `single_relations`，记录已完成事件表 `completed_events(event_index, finished_at, candidate_count)`。
   - **断点续跑机制**：
     - 若 `regulatory_network.duckdb` 已存在，先比对 `network_metadata` 中的配置指纹、输入指纹、代码指纹与当前运行一致；一致则读取 `completed_events` 跳过已完成事件继续运行，不一致则报错拒绝混写（需换新输出目录）。
     - 所有目标事件完成后统一执行 `rank_relations`，更新 `status = 'selected'` 并执行 `CHECKPOINT`。
   - **紧凑存储与导出**：
     - `single_relations` 不再重复存完整 `candidate_json` 与完整设计矩阵规格，仅存紧凑候选摘要（`route`、`priority`、`anchor_count`、`path_count`）与紧凑拟合系数/状态 JSON；`regulatory_network.tsv.gz` 只导出结构化数值与状态列（不含大段 JSON）。

---

### 2.4 核心建模、组成级联与支持域改造（`src/splicefate/model.py`、`src/splicefate/composition.py`）

1. **状态基清理与多上下文时间标准化**
   - `state_basis` 自动检测并移除与已有列（含截距）线性相关的重复指示列（例如 `stage` / `batch` / `branch` 共线时只保留一组独立列，记录实际保留列名于 `spec['state']`）。
   - 若数据包含多个独立轨迹上下文（如多阶段发育数据），伪时间 `z` 在各 `branch` 内部各自标准化，避免跨孤立轨迹混拼时间尺度。

2. **设计基改在观测细胞 `good` 上构造**
   - 在 `fit_design` 中，因子 `mean`/`scale`、样条节点 `knots`、样条中心化均值 `basis_mean` 以及状态正交投影 `state_projection` 均改在目标事件观测细胞 `good = (C.sum(1) > 0)` 上计算。
   - 交互项正交投影与 `composition` 投影保持在观测交集上计算；若某因子在 `good` 细胞内峰峰值 `< 1e-6`（无有效变异），标记 `kind = 'constant'` 并记录 `weak_identifiability = True`，避免构造空样条区间。
   - 单个单关系拟合若其唯一因子在 `good` 上为常数，直接标记 `constant_input_on_observed_cells` 跳过无效优化。

3. **固定离散度与分支数上限**
   - 保持 `model.dm_concentration: 20.0` 不变。
   - 新增配置项 `model.max_event_branches: 8`（整数 `>= 2`）。
   - 超过该分支数上限的事件在 `fit` 阶段直接记为状态 `too_many_branches_unmodeled`，写入 `all_event_status.tsv` 与 `fit_diagnostics.tsv`（GSE295932 全量中仅涉及 18 个超复杂事件）。

4. **组成级联低计数平滑（`composition.py`）**
   - `next_inputs` 更新上游 RBP 自身组成时，改用已配置的 `count_shrinkage`（默认 `2.0`）和该事件的全局经验比例 `pooled` 构造平滑基线：
     `P_ref = (C + s * pooled) / (N + s)`
   - 再在 logit 空间叠加相对 WT 的预测变化：
     `updated = softmax(log(max(P_ref, 1e-12)) + delta)`
   - 仍只更新有来源事件观测的细胞（`observed`）；当 `delta == 0`（含 `q = 0`）时严格短路返回原始输入，保证 `q = 0` 精确恒等。

5. **支持域修复、三项贡献拆分与悬空配置清理**
   - **联合支持域（`fit_event` / `predict_change`）**：
     - 真正启用 `model.support_neighbors`（默认 `10`，当观测细胞数不足时取 `min(k, n_good - 1)`）计算第 k 近邻距离的 `support_distance_quantile` 分位数作为局部密度支持半径。
     - 支持空间由“仅表达”扩展为 **`(表达 + 组成特征 + 标准化时间)`** 联合向量（仅包含本模型实际使用的特征，全部在 `good` 细胞上标准化）。
     - `input_outside_WT_range` 改用 `good` 观测细胞上的 WT 范围；额外增加 `input_below_observed_WT_p05`（或低支撑稀疏区标记）用于提示 KO 落入 `good` 细胞极低表达尾区的情况，不改变预测值，不做硬性幅度截断。
   - **三项贡献拆分**：
     - 在模拟阶段对每个受影响事件记录并输出：
       - `delta_expression_only`（第 1 轮仅施加直接表达剂量变化、冻结组成项的效应）
       - `delta_direct_composition`（第 1 轮同时施加直接目标表达与直接目标组成缩放相对仅表达项的增量）
       - `delta_cascade`（最终轮次相对第 1 轮直接效应的多轮传播增量）
       - 满足 `delta_probability == delta_expression_only + delta_direct_composition + delta_cascade`。
   - **悬空配置处理**：
     - 启用 `simulation.divergence_limit`：在多轮级联中，若相邻轮次最大 logit 差或状态更新范数超过该上限，立即熔断停止并记录 `stop_reason = 'divergence_limit_exceeded'`，避免发散。
     - 从配置与文档中清理未使用的 `simulation.gene_limit`。

---

### 2.5 模型、结构与模拟输出存储重构（`src/splicefate/cli.py`、`src/splicefate/simulation.py`、`src/splicefate/structure.py`）

1. **模型与结构去大 JSON**
   - `models` 产物拆为：
     - `models.duckdb`（或紧凑元数据表）：每个事件的 `design`、`coefficients`、`candidates` 摘要、`input_diagnostics`、状态与统计字段（不再嵌全细胞大数组与冗余 `anchors`/`paths`）。
     - `model_arrays.npz`：按紧凑连续数组与偏移量存储全部已拟合事件的 `WT_prediction`（`float64`）、支持域参考矩阵与索引；`observed_denominators` 直接从输入 `branch_counts.npz` 按需读取，不再冗余复制。
   - `structures` 产物拆为：
     - `structures_meta.json`：仅保存每个基因的 `chain_ids`、`event_branch_maps`、`retained_fraction`、`scope` 及列偏移量（不含计数矩阵）。
     - `structure_counts.npz`：所有建模基因保留剪接链合并为单个稀疏 `CSR` 矩阵（`cells × total_retained_chains`）。
   - `fit_summary.json` 校验新产物的 SHA256。

2. **流式模拟与按基因并行结构重建**
   - 新增配置项 `simulation.write_unchanged_cells: false`。
   - `simulate` 取消在内存中构造千万级 `list[dict]`，改为：
     - 逐事件（或分块）计算细胞级概率与三项贡献拆分，边算边累加 `event_stage_summary` 统计量；
     - 细胞级 `events.tsv.gz` 默认仅流式写出受影响事件（`target_relation_status` 为 `modeled_direct` / `modeled_cascade_only` 或 `max|Δp| > 0`），当 `write_unchanged_cells: true` 时流式写全量；
     - `structure_predictions` 仅对发生事件变化的基因按基因多进程并行执行 `reconstruct`，结果流式写入 `chains.tsv.gz` 和 `structure_diagnostics.tsv`。

3. **输出状态区分与三种评测聚合口径**
   - `event_stage_summary.tsv` 始终覆盖**全部事件 × 阶段 × 分支**，新增字段：
     - `target_relation_status`：`modeled_direct` / `modeled_cascade_only` / `not_selected_by_quality` / `no_prior_path` / `event_unmodeled`
     - `network_best_rank`、`network_best_selection_score`、`network_best_oob_gain`（当目标在单关系网络中存在但被筛除时可直接溯源）
     - 三种聚合口径：
       - `delta_probability_observed_cells`：仅在该事件 WT 有分子观测的细胞上均值
       - `delta_probability_molecule_weighted`：按 WT 事件观测分子数 `N_i` 加权的 pseudobulk 变化（`sum(N_i * p_pert_i)/sum(N_i) - sum(N_i * p_wt_i)/sum(N_i)`），作为与真实 KO 对比的主口径
       - `delta_probability_all_cells`：全细胞无权条件均值（保留兼容，原 `delta_probability` 列默认对齐主口径或显式标注并在 `evaluate` 中可指定列名）
     - 三项贡献汇总：`delta_expression_only_mw`、`delta_direct_composition_mw`、`delta_cascade_mw`

---

### 2.6 真实 KO 外部评测脚本与健壮性修复

1. **GSE295932 真实 KO 同坐标投影脚本（新增 `scripts/project_gse295932_ko_truth.py`）**
   - 从 `test/GSE295932/files/` 中读取与主数据集同条件（`GSM8961554` 刺激态）的非 Nontarget CRISPR KO 细胞（按官方 barcode 注释分组）。
   - 使用与 WT 完全相同的 `parse_blocks`、`groups.json` 公共锚点与分支定义，逐 KO 靶点统计每个事件分支在 WT 与 KO 细胞中的真实分子计数与 pseudobulk 比例变化 `delta_probability_truth`（附带 WT/KO 观测细胞数、总分子数与过滤审计）。
   - 单独对 `BUB3`（不在 1 跳剪接因子网内）与 `VIM`（阴性对照）输出外部覆盖与阴性对照审计表。
   - 输出标准 TSV，直接输入 `splicefate evaluate`（`evaluate` 支持 `--metric` 选择 `delta_probability_molecule_weighted` / `delta_probability_observed_cells` / `delta_probability_all_cells`，默认 `delta_probability_molecule_weighted`），全程只做外部评测，不参与任何训练或选参。

2. **健壮性修复**
   - **全量 `execution.bootstrap` 修复（`src/splicefate/cli.py`）**：改为逐事件收集收敛成功的 bootstrap 样本，不再因全部 2 万余事件中任意 1 个事件未收敛就丢弃整轮结果；每个事件在输出中明确记录各自的 `bootstrap_successful_fits`。
   - **`--check` 补齐 R 环境校验（`src/splicefate/runner.py`、`src/splicefate/cytotrace.py`）**：当 `trajectory.root_method == 'cytotrace'` 时，`--check` 实际验证 `Rscript` 可执行，且 `Matrix`、`HiClimR`、`ccaPP`、`nnls` 四个 R 包均可成功 `requireNamespace`，并校验官方 `CytoTRACE.R` SHA256。

---

## 3. 新主配置方案（`configs/gse295932_stimulated_hg38.yaml`）

保留原 `configs/gse295932_nontarget_hg38.yaml` 不覆盖，新增拆开后的刺激态主配置 `configs/gse295932_stimulated_hg38.yaml`（并同步更新通用字段到其余配置）：

```yaml
project_root: ..
input:
  genes: data/processed/GSE295932_Nontarget_stimulated_hg38_cluster/genes
  events: data/processed/GSE295932_Nontarget_stimulated_hg38_cluster/events
  database: data/external/human_rbp_unified/human_rbp.duckdb
  fasta: ../reference_documents/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
  velocity_times: null
  intron_counts: data/processed/GSE295932_Nontarget_stimulated_hg38_cluster/regions
assembly: hg38
output: runs/GSE295932_stimulated_v1_full
runtime:
  python_binary: python
  threads: 24
execution:
  stages: [prepare, fit, simulate]
  event_limit: null
  bootstrap: 0
clustering:
  enabled: false
trajectory:
  S_calibration: intron_phase
  method: dpt
  root_method: cytotrace
  isolate_stages: [stimulated]
priors:
  targets: [ALL]
  partner_hops: 1
  extra_targets: []
  include_family_shared: true
  max_partners_per_factor: 8
  min_experiment_score: 400
  allow_unspecified_association: false
  window_bp: 100
  motif_pvalue: 0.0001
  min_event_cells: 8
  min_event_molecules: 20
model:
  response: spline
  spline_knots: 3
  state_ridge: 1.0
  ridge: 10.0
  interaction_ridge: 30.0
  dm_concentration: 20.0
  max_event_branches: 8
  max_factors: 5
  network_selection:
    method: all_single_relations_then_weighted_quality
    min_response_score: 0.0
    min_oob_fit_gain: 0.0
    min_positive_gain_fraction: 0.6
    min_stability: 0.3
    bootstrap_repeats: 10
    bootstrap_seed: 42
    weights:
      fit_gain: 0.45
      stability: 0.45
      response: 0.10
  support_neighbors: 10
  support_distance_quantile: 0.95
  count_unit: integer_umi
simulation:
  target_sets:
    - [CELF2]
    - [CWC27]
    - [HNRNPD]
    - [HNRNPLL]
    - [PCBP2]
    - [RBFOX2]
    - [SFPQ]
    - [SNRPD2]
    - [SRRM1]
    - [SRSF5]
    - [SRSF9]
  strengths: [0.0, 0.5, 1.0]
  max_rounds: 3
  convergence_tolerance: 1.0e-06
  divergence_limit: 10.0
  write_unchanged_cells: false
  max_chain_states: 128
  structure_lambda: 10.0
  structure_prior_strength: 1.0
```

---

## 4. 实施顺序与验证计划

1. **阶段 A：数据准备与外部真值脚本**
   - 升级 `scripts/convert_gse295932.py` 与 `scripts/quantify_gse295932_regions.py` 支持 `--condition stimulated`，生成 `GSE295932_Nontarget_stimulated_hg38_cluster/` 及其 `regions/`。
   - 新增 `scripts/project_gse295932_ko_truth.py`，生成与刺激态 `groups.json` 严格对齐的外部 KO 评测表（不触碰训练输入）。

2. **阶段 B：先验、存储与核心模型改造**
   - `io.py`：更新配置校验（`partner_hops`、`extra_targets`、`max_event_branches`、OOB 门槛、`write_unchanged_cells`，移除 `gene_limit`）。
   - `priors.py`：实现 `partner_hops` 解耦与 `priors_meta.json` + `priors.duckdb` 存储。
   - `model.py` & `composition.py`：实现共线列清理、`good` 观测细胞设计基、分支数上限、低计数级联平滑、联合支持域（`support_neighbors`）与三项贡献拆分接口。

3. **阶段 C：网络并行、OOB 评分与流式模拟**
   - `network.py`：实现 OOB 增益与稳定性评分、事件级基线缓存、24 核按事件多进程并行、DuckDB 按事件断点续跑与紧凑导出。
   - `structure.py`、`simulation.py`、`cli.py`、`runner.py`：实现 `models.duckdb` + `model_arrays.npz`、`structure_counts.npz`、流式模拟写盘、按基因并行结构重建、5 类 `target_relation_status`、3 种聚合口径、`execution.bootstrap` 修复及 `--check` R 环境检查。

4. **阶段 D：测试、小规模真实验证与文档同步**
   - 补齐本地单元测试（覆盖 OOB 无自由度偏置、事件级基线缓存数值恒等、断点续跑一致性、DuckDB/NPZ 读写、`good` 设计基、低计数平滑、三项贡献加和恒等、联合支持域、5 类状态码、3 种聚合口径、R 依赖 `--check`）。
   - 使用刺激态完整 WT（662 细胞）+ 少量事件（如 20 事件）跑通完整 `prepare → fit（含断点续跑演练）→ simulate（11 KO × 3 剂量）→ evaluate`，验证耗时、内存占用、磁盘体积及全部数值不变式。
   - 同步更新 `README.md`、`docs/IMPLEMENTATION_v0.3.md`、`docs/REGULATORY_NETWORK_LEARNING_20261121.md`，核对无本地路径与测试文件后提交代码。
   - 经用户确认后，再启动 `configs/gse295932_stimulated_hg38.yaml` 的第一轮完整全量运行。
