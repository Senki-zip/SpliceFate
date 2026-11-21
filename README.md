# SpliceFate

SpliceFate 使用 WT 单细胞表达和长读长剪接计数，结合 motif、CLIP、蛋白互作和复合体信息，预测 RBP 扰动后的剪接选择与剪接链组成变化。

核心模型是状态条件下的样条响应和 Dirichlet–multinomial 复合计数目标。RBP 总表达与自身局部剪接组成共同作为输入；组成作用系数从 WT 学习，并用于同步多轮传播。

## 安装

```bash
conda env create -f environment.yml
conda activate SpliceFate
python -m pip install -e '.[conversion]'
```

CytoTRACE 起点选择使用原版 0.3.1 R 源码，位于 `tools/CytoTRACE_official/`。R 运行依赖为 R≥4.0、Matrix、HiClimR、ccaPP、nnls。

```r
install.packages(c("Matrix", "HiClimR", "ccaPP", "nnls"))
```

Python 和 R 默认从当前环境的 PATH 查找。跨环境运行可设置 `SPLICEFATE_PYTHON`、`SPLICEFATE_RSCRIPT`，或修改 YAML 中的 `runtime.python_binary`、`trajectory.cytotrace.rscript_binary`。备份默认写入项目 `backups/`，可通过 `SPLICEFATE_BACKUP_DIR` 指定目录。

## 配置和启动

```bash
# 默认使用 configs/gse295932_nontarget_hg38.yaml
python scripts/run_splicefate.py --check
python scripts/run_splicefate.py

# 选择其他数据配置
python scripts/run_splicefate.py --config configs/nemo_p21_v0.3.yaml
```

所有运行参数写入 YAML。相对路径以配置中的 `project_root` 为基准。

```yaml
project_root: ..
input:
  genes: data/processed/GSE295932_Nontarget_hg38_cluster_20261121/genes
  events: data/processed/GSE295932_Nontarget_hg38_cluster_20261121/events
  database: data/external/human_rbp_unified/human_rbp.duckdb
  fasta: ../reference_documents/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
  velocity_times: null
  intron_counts: null
assembly: hg38
output: runs/GSE295932_Nontarget_cluster_RBFOX2_20261121
runtime:
  python_binary: python
  threads: 24
execution:
  stages: [prepare, fit, simulate]
  event_limit: null
  bootstrap: 0
clustering:
  enabled: false
priors:
  targets: [ALL]
simulation:
  target_sets: [[RBFOX2]]
  strengths: [0.0, 0.5, 1.0]
```

以上展示主要运行字段，完整配置见 [人 Jurkat 配置](configs/gse295932_nontarget_hg38.yaml)和[小鼠配置](configs/nemo_p21_v0.3.yaml)。

- `execution.stages`：按顺序执行 prepare、fit、simulate。复用已拟合结果时设为 `[simulate]`。
- `execution.event_limit`：`null` 拟合全部合格候选事件；正整数限制拟合事件数。
- `execution.bootstrap`：WT 细胞在 stage×batch 内重采样，重拟合并运行级联的次数。
- `runtime.threads`：BLAS、OpenMP、Numba 和 Scanpy 的线程预算。CytoTRACE R 进程数由 `ncores` 单独配置。
- `priors.targets`：候选网络的起始因子；`[ALL]` 选择全部已表达剪接相关候选。
- `simulation.target_sets`：每个内层列表定义一次干预；`[[DDX5], [DDX5, RBFOX2]]` 分别表示单独和联合扰动。
- `simulation.strengths`：0=无扰动、0.5=降低一半、1=完全 KO。脚本执行靶点组合与强度的笛卡尔积。

prepare 和 fit 各执行一次，随后依次执行模拟。输出保存配置快照、输入和代码指纹、线程信息、运行计划及阶段状态。改变输入、模型设置或代码时使用新的输出目录。

原有分阶段 CLI 继续可用：

```bash
splicefate prepare --config configs/gse295932_nontarget_hg38.yaml
splicefate fit --config configs/gse295932_nontarget_hg38.yaml
splicefate simulate --config configs/gse295932_nontarget_hg38.yaml --targets RBFOX2 --strength 1
```

## 输入格式

矩阵均为 **cell×feature**，以 TSV 中的行序定义轴。默认计数为非负整数。

| 目录 | 文件 | 内容 |
|---|---|---|
| genes | `cells.tsv` | 唯一 cell ID，以及 stage、stage_order、batch 等元数据 |
| genes | `genes.tsv` | 稳定 gene ID、gene_name，及可选坐标注释 |
| genes | `gene_counts.npz` | 完整基因整数计数，用于表达输入和文库分母 |
| genes | `isoforms.tsv` | isoform、gene、JSON structure 与 exons |
| genes | `isoform_counts.npz` | 实际观测剪接链的计数 |
| genes | `quantification.json` | 计数单位、来源、assembly 和验收信息 |
| events | `cells.tsv` | 与 genes 相同的细胞集合 |
| events | `genes.tsv` | Seurat v3 选择的 3,000 HVG |
| events | `groups.json` | 事件 ID、gene、chromosome、strand、公共 anchors 和各分支 junctions |
| events | `branch_counts.npz` | 事件各分支的计数 |
| events | `branch_offsets.npy` | 每个事件分支的列区间 |
| events | `normalization_factors.tsv` | 完整文库大小及 CP100000 比例 |

```json
{"isoform":"SC_example","gene":"ENSG00000000001","structure":["junction:chr1:100-200:+"],"exons":[[50,100],[200,300]]}
```

坐标统一为 **0-based half-open**，支持 mm10 和 hg38。motif 扫描、CLIP 与事件使用相同的原生 assembly。细胞 ID 建议采用 `sample:barcode`。

事件由公共锚点覆盖下的竞争分支定义，N=0 表示该事件未观测。分支顺序由 `groups.json` 的结构决定。多个事件可共享来源分子，拟合目标采用复合计数形式。剪接链保存实际观测结构及终端覆盖范围。

## GSE295932 数据转换

```bash
python scripts/download_gse295932.py
python scripts/separate_gse295932_controls.py
python scripts/convert_gse295932.py
python scripts/build_human_model_priors.py
```

下载器使用 `test/GSE295932/metadata/` 中的 GEO 文件清单与样本元数据。数据来源见 [GEO GSE295932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295932)及[论文](https://doi.org/10.1038/s41467-025-60902-2)。

当前训练输入为全转录组 Nontarget 对照：未刺激 387、刺激后 662，共 1,049 个细胞。刺激条件保留独立轨迹上下文。

| 字段 | 单位 | 当前选择 |
|---|---|---|
| `size` | 支持 reads，包含 PCR 扩增 | 来源审计 |
| `cluster` | 校正结构对应的整数 UMI 簇，散射／小簇过滤前 | **G、链、事件共同使用** |
| `count` | 过滤后的加权计数，可为分数 | 历史版本比较 |

当前整数版有 16,699 个基因、217,053 条剪接链、28,266 个事件和 60,040 个分支。G 总数 5,686,581，HVG 按当前 G 重新计算。cluster 保留过滤前的小簇，其噪声差异在与加权版的比较中记录。

输出目录为 `data/processed/GSE295932_Nontarget_hg38_cluster_20261121/`。来源、字段选择、坐标核验和两版比较见 [开发记录](docs/GSE295932_COUNT_SELECTION_20261121.md)。历史分数模式可通过 `model.count_unit: author_fractional_umi_estimates` 和配套 quantification 记录启用，对应 gamma 扩展 DM 目标。

## 小鼠 AllInfo＋MEX 转换

```bash
splicefate-convert --input test/processed --output data/processed/NEMO_P21_author_20261117 --gtf ../reference_documents/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz
```

小鼠 G 来自配对短读长整数 UMI，S 来自分子级 AllInfo。按 sample、barcode、gene、exact UMI 去重，冲突结构单列审计。局部竞争包含 cassette、共享 donor 和共享 acceptor，两侧公共锚点各 10 bp。事件保留门槛为总分子≥10、有效细胞≥3、至少两个分支被观测。

NEMO P21 转换产物包括 35,317 个细胞、31,053 个基因、3,000 HVG、147,447 个事件和 359,082 个分支。

## 模型流程

1. **表达参考**：完整 G 文库归一化到 CP100000/log1p，使用固定 Seurat v3 HVG 建立 PCA 和邻接图。Leiden 为可选步骤；关闭时使用 global 上下文，保留显式条件和连通分量边界。
2. **CytoTRACE→DPT**：每个上下文使用完整基因计数运行官方 CytoTRACE 0.3.1，最高有效分细胞作为起点；已知多阶段时在最早阶段内选择。评分同分按 cell ID 排序。
3. **S 局部校准**：DPT 为初始轴，以可靠 E/I 相位或剪接组成约束做局部小幅修正。全部事件使用同一份时间；缺少 S 校准输入时使用原 DPT。
4. **候选关系**：从同物种 motif、CLIP、PPI 与复合体生成直接及伙伴介导的事件候选，并保存证据来源。
5. **先学习网络、再筛选输入**：对全部有先验证据的候选 RBP–事件关系逐条拟合条件 DM 响应，包含 RBP 总表达及可用的自身局部组成。完整单关系网络保存后，按 WT 拟合响应变化评分，为每个事件选择最多 5 个因子，再联合拟合阶段、批次、轨迹背景、总表达、交互与组成项。
6. **虚拟扰动**：在表达原单位上乘以 `(1-q)`，转换回 log1p 后计算相对 WT 的响应差分。联合干预同时求值。
7. **组成级联**：预测改变 RBP 自身剪接状态后，同步更新其组成输入，再由冻结的响应模型计算下一轮。
8. **结构输出**：将最终事件响应约束到实际观测链的概率 simplex，记录组成、保留结构比例及优化状态。

输出重点为事件分支概率变化、剪接链使用比例变化和 WT 学习的组成响应系数。独立扰动评价通过匹配目标、剂量、事件、条件和分支进行。

```yaml
model:
  max_factors: 5
  network_selection:
    method: all_single_relations_then_top_response
    min_response_score: 0.0
```

`max_factors` 限制最终联合模型的输入数。完整单关系学习覆盖所有候选，按 CausalBridge 的“单关系学习→WT 响应评分→组合拟合”顺序执行。评分使用同一上下文、阶段和批次内相邻 WT 细胞的因子响应变化；状态背景保持配对，事件缺失观测保留原位置。完整定义及保存格式见[网络学习说明](docs/REGULATORY_NETWORK_LEARNING_20261121.md)。

## 异构体组成和多轮传播

RBP 自身的公共锚点事件定义局部异构体分组。默认每因子选择一个计数有支持的事件，最多四分支；按少数分支分子支持排序。

```yaml
model:
  isoform_state:
    enabled: true
    max_events_per_factor: 1
    max_branches: 4
    min_branch_molecules: 20
    min_branch_cells: 10
    min_joint_cells: 20
    count_shrinkage: 2.0
    ridge: 10.0
```

组成输入为 `x × (p_k-pooled_k) × N/(N+count_shrinkage)`，其中 x 是 log1p 归一化总表达，p 是实际观测分支比例。K 个分支使用 K−1 个对比，总表达项独立保留。每条下游关系学习有符号系数，可包含时间交互。

组成项按实际观测门控：未观测到上游状态的细胞使用总表达响应；联合支持用于确定组成项的可拟合性。因子自身事件的预测使用其他输入，避免将目标分支直接作为自身预测变量。这个规则只作用于响应模型，伪时间使用共享全部特征。

每轮相对同一 WT 基线重新求值，直接靶点持续保持指定剂量。组成更新仅用于有来源观测的细胞。达到轮数预算或输入变化阈值后停止；轮数表示传播层次。开启 bootstrap 时，对 WT 重采样重拟合并重新执行整条级联。

## S 校准输入

`trajectory.S_calibration` 支持 `intron_phase`、`composition` 和 `hybrid`。内含子模式需要独立的 cell×gene `exon_counts.npz`、`intron_counts.npz`、细胞完整文库分母、基因轴和 quantification 元数据。

先估计基因参考关系 `I ≈ slope × E`，相位残差与邻居外显子状态差提供局部方向，随后限制 DPT 修正幅度、固定起点和阶段边界。没有支持的细胞保持原值。

当前 NEMO 和 GSE295932 输入以剪接链／事件计数为主，独立 E/I 矩阵尚未量化。AllInfo 的 intron chain 是已切除的缺口；比对块中的 intronic coverage 可以用于后续区域定量。`intron_counts: null` 时使用 CytoTRACE 定向的 G-only DPT。

可选完整动力学比较接口 `splicefate.velocity` 使用经分类核验的 gene-level spliced/unspliced 整数层，输出 `cell`、`latent_time`、`qc_pass`；路径写入 `input.velocity_times`。

## RBP 数据库

人源模型先验位于 `data/external/human_rbp_unified/`：native human CISBP/ATtRACT/RBP2GO、hg38 ENCORI CLIP、STRING physical 和 CORUM/ComplexPortal。原始来源与构建日期保留在数据库元数据和清单中。

小鼠主库位于 `data/external/mouse_rbp_unified_20261116/`，其外部位点文件位于同目录 `sites/`。蛋白维度扩充包括 UniProt 序列、功能区域、InterPro/Pfam 坐标、RBD 注释、结合肽段、参考异构体及经翻译核对的 domain—外显子映射。

```bash
python scripts/query_mouse_rbp_database.py --factor Ptbp1 --protein-only
python scripts/collect_mouse_protein_dimensions.py
python scripts/integrate_mouse_protein_dimensions.py
python scripts/finish_mouse_protein_dimensions.py
```

蛋白资源下载使用项目内 REST 客户端。结构域、家族和超家族类型分别保存；不同来源的坐标匹配保留来源与序列核验结果。家族共享 motif 采用同物种 InterPro Type=Family 一跳扩展，标注 `family_shared_candidate`。蛋白注释用于结果解释，组成级联直接从 WT 计数学习权重。

## 输出

| 文件 | 内容 |
|---|---|
| `run_status.json`、`pipeline.log` | 阶段、任务进度与日志 |
| `launch_config.resolved.yaml`、`runtime.json` | 运行参数与线程信息 |
| `provenance.json` | 输入、代码和准备产物校验值 |
| `cells.tsv`、`trajectory_diagnostics.json` | 细胞元数据、CytoTRACE、起点与 DPT/GS 时间 |
| `prepared_arrays.npz` | 全事件共享的一维时间、PCA 和因子输入 |
| `priors.json` | 因子轴、锚点、伙伴路径和局部组成状态 |
| `models.json`、`fit_summary.json` | 拟合模型与统计 |
| `regulatory_network.duckdb` | 全部单关系的拟合曲线、证据、得分及筛选结果 |
| `regulatory_network.tsv.gz`、`regulatory_network.summary.json` | 完整关系表与网络统计 |
| `isoform_response_weights.tsv` | 上游组成对下游事件的响应系数及联合支持 |
| `all_event_status.tsv`、`fit_diagnostics.tsv` | 事件计数、拟合与支持诊断 |
| `simulations/*/events.tsv.gz` | 各细胞分支概率及变化 |
| `simulations/*/chains.tsv.gz` | 各细胞剪接链比例及变化 |
| `simulations/*/structure_diagnostics.tsv` | 结构约束、冲突和优化状态 |
| `simulations/*/simulation.json` | 传播历史和停止原因 |

## 开发记录和参考

- [GSE295932 计数选择](docs/GSE295932_COUNT_SELECTION_20261121.md)
- [模型实现](docs/IMPLEMENTATION_v0.3.md)
- [伪时间与级联设计](docs/SPLICEFATE_PSEUDOTIME_AND_CASCADE_v0.3.md)
- [CausalBridge](https://github.com/Senki-zip/CausalBridge_en)：候选路径与冻结 WT 差分参考。
- [CytoTRACE](https://doi.org/10.1126/science.aax0249)：DPT 起点评分。
- [scVelo](https://github.com/theislab/scvelo) 和 [Region Velocity](https://doi.org/10.1093/nar/gkad969)：RNA 相位方向与动力学参考。

第三方软件条款见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
