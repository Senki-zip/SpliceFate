# SpliceFate v0.3

WT 条件下的剪接扰动建模：G 参考聚类、S 辅助伪时间、motif/CLIP 与蛋白互作先验、计数响应、合法剪接链组成和WT 学习的剪接组成级联。

当前为可执行研究实现，尚未经过独立 KO/KD 效果验证。当前不输出已校准因果效应、成熟 RNA 总量变化或命运概率。实现与特定数据集分开；没有要求固定 161 个细胞或指定阶段数量。

## 小鼠蛋白／结构域扩充库

2026-11-20 开始在原 `data/external/mouse_rbp_unified_20261116/mouse_rbp.duckdb` 增加蛋白层；原 motif、CLIP、家族共享与互作表不重建。目录日期保留，实际来源日期、序列版本、文件 SHA256 另行记录。新增原始资源、下载状态、验收及覆盖表位于 `data/external/rbp_protein_dimensions_20261120/`。

范围为有来源支持的小鼠 RBP，及 GO 剪接相关候选（包括不直接结合 RNA 的因子）。UniProtKB 提供完整蛋白记录、功能区域、活性／结合位点、定位、可变剪接和蛋白异构体 FASTA；InterPro API 提供 integrated 与 Pfam 等成员库的坐标匹配；RBP2GO 已有的结构域、无序区比例和实验结合肽段被规范化；本地 mm10/mm39 GTF 与 FASTA 用于编码区核验。

```bash
# 仅查蛋白／domain／异构体／编码区，不扫描整个 motif 位点库。
python scripts/query_mouse_rbp_database.py --factor Ptbp1 --protein-only --limit 20

# 可恢复下载；已验证文件复用。UniProt 通过 life-science-research 技能客户端获取。
python scripts/collect_mouse_protein_dimensions.py

# 完整来源下载后整合；--allow-partial 明确表示增量快照，不能称为全来源完成。
python scripts/integrate_mouse_protein_dimensions.py

# 后台下载期间定期入库；下载全部完成后再完整验收并备份。
python scripts/finish_mouse_protein_dimensions.py
```

新增主要查询层：

| 表／视图 | 内容 |
|---|---|
| `protein_entries`、`protein_sequences` | 小鼠蛋白身份、序列、版本和序列校验值 |
| `protein_factor_mapping`、`protein_accession_resolution` | 原库因子关联、主／次 accession 解析、歧义与基因冲突标签 |
| `domain_catalog`、`factor_domain_annotations` | 结构域词典、原 RBP2GO 的 domain/RBD 数量与归属；没有坐标的记录不伪造位置 |
| `protein_domain_matches`、`protein_domain_fragments` | UniProt、InterPro/Pfam 等来源的位置，保留不连续片段；Family/Superfamily/Domain 分开 |
| `sequence_verified_domain_matches`、`factor_protein_domains` | 仅序列核对一致的 domain/repeat/zinc-finger 坐标和解析明确的因子关联 |
| `protein_features`、`protein_feature_evidence`、`protein_function_comments` | 定位、功能区域、结合／活性位点、PTM、变异／突变注释与 ECO/文献来源 |
| `protein_isoforms`、`protein_isoform_edits`、`protein_isoform_feature_consequences` | 参考蛋白异构体及经完整序列核验的区域保留、删除、重叠改动；不是功能活性分数 |
| `rna_binding_peptides`、`protein_peptide_domain_overlaps` | 实验结合肽段、原蛋白坐标、序列一致性和 domain 重叠；重叠不等于 domain 特异结合验证 |
| `reference_coding_segments`、`reference_coding_transcripts` | 本地两个参考版本的 CDS/UTR/起止密码子、相位及翻译完整性 |
| `reference_transcript_function_context` | 参考转录本类型及来源已有的 NMD 标签；不把标签自动转移到新观测链 |
| `protein_transcript_links`、`protein_transcript_link_audit` | 来源转录本／蛋白交叉引用、版本和本地参考的匹配状态 |
| `protein_local_transcript_matches`、`protein_domain_genomic_blocks` | 同基因且 CDS 翻译与蛋白序列完全一致的转录本关联；domain 映射回有参考支持的外显子区间 |
| `protein_dimension_coverage`、`protein_dimension_policy`、`protein_dimension_input_assets` | 因子覆盖、解释边界、来源与校验信息 |

所有新增坐标为 **0-based half-open**；来源 1-based 坐标保留在原始文件／JSON 中。不同数据库在同一位置的多个匹配是重复来源注释，不能直接相加为生物学 domain 数。InterPro 序列与 UniProt 不一致时保留审计，坐标不跨序列转移；蛋白异构体只在来源剪接编辑能够完整重建其 FASTA 时投射功能区域。基因名相似不会触发新蛋白映射，歧义 accession 和跨基因冲突不进入已解析视图。

当前扩充没有给 domain 缺失自动赋予 0–1 功能权重，也没有把 mutagenesis、疾病变异或蛋白功能描述导入 WT 扰动响应先验。参考转录本映射不等于补全本数据的截短剪接链；新观测链仍需独立判定 CDS 与覆盖完整性。

`integration_status.json` / `validation.json` 是新增层的验收；`download_status.json` 是外部来源下载状态；`pipeline_status.json` 是持续整合流程状态。部分来源仍在下载时，已验证快照可查询，但必须保留 `partial_sources_integrated` 标签。数据库文件更新后，已有模型的输入指纹会变化，后续建模应在新输出目录重新 prepare。

## 安装与命令

推荐使用配置启动入口：

```bash
# 默认读取 configs/gse295932_nontarget_hg38.yaml，自动使用其中的 Python 解释器。
python scripts/run_splicefate.py

# 仅读取、检查配置和输入并打印任务计划，不执行模型或写入输出。
python scripts/run_splicefate.py --check

# 使用其他配置文件。
python scripts/run_splicefate.py --config configs/your_run.yaml
```

输入目录、数据库、FASTA 在 `input`，输出路径在 `output`；相对路径以配置中的 `project_root` 为基准，不受启动所在目录影响。所有运行选择写入 YAML：

```yaml
runtime:
  python_binary: /home/huangtao/miniconda3/envs/SpliceFate/bin/python
  threads: 4
execution:
  stages: [prepare, fit, simulate]
  event_limit: null
  bootstrap: 0
simulation:
  target_sets:
    - [Ptbp1]          # 单独扰动 Ptbp1
    - [Ddx5]           # 单独扰动 Ddx5
    - [Ddx5, Rbfox2]   # 同时扰动两个靶点
  strengths: [0.0, 0.5, 1.0]
```

启动脚本只接受配置选择和 `--check`，路径、线程、事件范围、bootstrap、KO 靶点及强度均从配置读取。`event_limit: null` 拟合全部满足候选与计数条件的事件；正整数用于限制事件数的检查。准备和拟合各执行一次，随后按 `target_sets × strengths` 逐项模拟；0=无扰动、0.5=降低一半、1=完全 KO。当前默认人源配置为 `priors.targets: [ALL]`、`clustering.enabled: false`，仅模拟 RBFOX2 的三种强度，共 3 个模拟任务；线程数保留配置现有值。

`clustering.enabled: false` 关闭 Leiden 计算、已有 Leiden 标签复用和 PAGA 合并，所有细胞分配到 `global` 上下文。PCA、邻接图、CytoTRACE 选根、DPT 和显式条件／连通分量边界保留；`leiden` 输出列写为 `all`。设置为 `true` 可恢复此前分群流程。

`configs/gse295932_nontarget_hg38.yaml` 是当前默认人 Jurkat 配置：RBFOX2、ALL 先验因子、关闭 Leiden，并显式分离 `unstimulated` / `stimulated` 条件。当前整数 cluster 转换产物在 `data/processed/GSE295932_Nontarget_hg38_cluster_20261121/`，人源 hg38 先验在 `data/external/human_rbp_unified/human_rbp.duckdb`；与原小鼠输入和数据库分开。小鼠运行显式指定 `--config configs/nemo_p21_v0.3.yaml`。

`priors.targets` 控制准备阶段的候选因子种子和伙伴范围，`simulation.target_sets` 控制实际要做的干预，二者用途不同。新增 KO 靶点时，应先确保它已包含在准备的候选因子范围；若要保证按该靶点生成候选路径，将它加入 `priors.targets` 并使用新输出目录。

`runtime.threads` 在导入数值库前设置 BLAS、OpenMP、Numba 等线程环境，也设置 Scanpy 的任务数；它控制线程预算，不表示同时启动多组事件拟合。官方 CytoTRACE 的 R 工作进程数仍由 `trajectory.cytotrace.ncores` 单独配置。配置的解释器至少需安装本项目依赖；最初启动脚本的 Python 只需可导入 PyYAML，随后自动切换至配置解释器。

输出目录保留 `pipeline.log`、`run_status.json`、`execution_plan.json`、`runtime.json` 和 `launch_config.resolved.yaml`。已有准备结果时，把 `execution.stages` 改为 `[fit, simulate]` 或 `[simulate]` 可以复用模型；不重复执行 prepare。修改线程、任务强度或已有候选范围内的 KO 组合不改变冻结模型指纹，输入、模型参数和代码仍执行一致性检查。新输入或新的建模配置使用新输出目录。原有 `splicefate prepare/fit/simulate/run` 命令也继续可用。

专用 Conda 环境：

```bash
conda env create -f environment.yml
conda activate SpliceFate
python -m pip check
python -m pytest -q
```

`environment.yml` 固定已验证的核心依赖版本，并包含 Seurat v3 转换依赖、PyArrow、pysam 与 Biopython。完整安装快照和检查结果保存在 `runs/environment_setup_20261119/`。默认 CytoTRACE 选根还需要下述独立 R 运行时；完整 scVelo 和 R/Region Velocity 不属于当前流程的必需依赖。

```bash
python -m pip install -e .

splicefate prepare --config configs/mouse_v0.3.yaml
splicefate fit --config configs/mouse_v0.3.yaml
splicefate simulate --config configs/mouse_v0.3.yaml --targets Ddx5 --strength 0.5
splicefate simulate --config configs/mouse_v0.3.yaml --targets Ddx5 Rbfox2 --strength 1

# 单次完成三个阶段；限制事件数时仅为实现检查，不是完整训练。
splicefate run --config configs/mouse_v0.3.yaml --output runs/example_smoke --event-limit 6 --targets Ddx5 --strength 0

# 可选 WT 细胞重采样区间；不使用外部 KO 标签拟合。
splicefate simulate --config configs/mouse_v0.3.yaml --targets Ddx5 --strength 0.5 --bootstrap 20

# 用固定外部真值评价，命令不训练、不选参。
splicefate evaluate --prediction event_stage_summary.tsv --truth heldout_perturbation.tsv --output benchmark.json
```

当前工作区专用解释器为 `/home/huangtao/miniconda3/envs/SpliceFate/bin/python`。核心测试为 `python -m pytest -q`，不依赖项目真实数据或联网。

真实数据安装试运行命令：

```bash
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 MPLBACKEND=Agg \
python scripts/run_model_smoke.py --config configs/nemo_p21_v0.3.yaml \
  --output runs/NEMO_P21_environment_smoke --events 12
```

该命令保留全部细胞，显式选取少量计数有支持的事件验证 prepare→fit→零剂量和部分剂量模拟，原转换数据和模型参数不变。它检查概率、链组成及严格 no-op，不代表完整训练或生物学 benchmark。先验 DuckDB、参考 FASTA 和配置路径必须存在；主库缺失时明确终止，不使用虚构先验替代。

`configs/mouse_v0.3.yaml` 是当前小鼠数据适配示例。更换数据时修改输入和输出路径，保持下述 schema；使用新输出目录。配置中的融合权重、正则、复杂度和传播预算是显式工程初值，尚未由 KO benchmark 优化。准备后改动代码、输入、配置或准备产物会触发一致性检查，不静默复用旧缓存。

## 输入契约

**当前文件格式以本节和转换命令为准。** 标准输入为 TSV/JSON 注释加 SciPy 稀疏 NPZ 矩阵；不能把 FASTQ、BAM、作者 AllInfo、TPM/PSI 或 log 表达直接放到计数矩阵路径。全部矩阵方向是 **行=细胞，列=特征**，当前默认计数为非负整数，GSE295932 使用作者 `cluster` 列；历史分数估计模式必须显式选择，不在当前配置启用；`.npy` 仅用于事件分支 offset。

### GSE295932 人 Jurkat Nontarget 转换

```bash
python scripts/convert_gse295932.py
# 默认使用本地 GRCh38-2024-A 的 GENCODE v44 注释。
python scripts/build_human_model_priors.py
python scripts/run_splicefate.py --check
```

转换只读分离好的全转录组 Nontarget：未刺激 387、刺激后 662，共 1,049 个细胞；不纳入 KO、靶向富集或单独的 ONT/PacBio WT 验证文库。复用 `splicefate.convert` 的剪接链标识、公共锚点事件发现与分支统计，适配 Longcell 的空格分隔、带引号表格。作者 exon intervals 为 1-based inclusive，转换到 0-based half-open；对 GENCODE v44 的抽查中，start−1 匹配 70,864 个 junction，未减一匹配 0 个。

**当前选择作者 `cluster` 整数列作为 G 和 S 的统计单位。** `size` 是含 PCR 重复的 reads 数；`cluster` 是分配到校正剪接结构的 UMI 簇数（小簇散射过滤之前）；`count` 是过滤后的加权估计，可能为分数。当前不使用 size 或 count，也不对小数取整。未刺激 cluster 总数 2,087,100，刺激后 3,599,481，共 5,686,581。`count=0` 但 cluster>0 的来源行按其 cluster 保留。选择细节、来源代码与局限见 [开发文档](docs/GSE295932_COUNT_SELECTION_20261121.md)。

G 用全转录组长读长的每细胞／基因 cluster 求和，不声称是配对短读长 G；分母为完整 gene cluster library。S 用同一计数列构建观测链和公共锚点事件分支，int64 保存，`model.count_unit: integer_umi` 对应整数 DM 复合计数目标。HVG 依本次 G 重新计算 Seurat v3 3,000 个。

当前输入目录 `data/processed/GSE295932_Nontarget_hg38_cluster_20261121/`，旧 `data/processed/GSE295932_Nontarget_hg38/` 作为加权 count 历史版保留，不再被默认配置读取。cluster 仍是作者处理后的 UMI 簇，且缺少作者散射过滤，可能保留错误小簇，不宣称是未经处理或逐分子验证的真实 UMI whitelist。分数估计兼容模式仍可显式选择并提供 quantification.json，但不是本次建模输入。

输出包含标准 `genes/`、`events/`、Seurat v3 3,000 HVG、CP100000 分母、来源 SHA256、`conversion_report.json`、`validation.json`、`model_input_check.json`。`conversion_staging.duckdb` 是可再生的中间表，不是模型必需文件。未知参考基因保留在 G 并单列审计；只有身份与坐标可解释的基因生成 S。

人源先验从已下载多物种 atlas 中提取 native human CISBP/ATtRACT/RBP2GO、hg38 ENCORI CLIP，并补充人 STRING physical 及 human CORUM/ComplexPortal；无小鼠投射。家族共享仍按同物种 InterPro Type=Family 一跳扩展并标记。当前没有导入 human IntAct，不将空缺当作已有证据。基因组参考已在本地，不重复下载。

### 从本次下载的作者文件转换

```bash
python -m pip install -e '.[conversion]'

python -m splicefate.convert \
  --input test/processed \
  --output data/processed/NEMO_P21_author_20261117 \
  --gtf /path/to/mm10/genes.gtf.gz
```

命令会流式读取四份 AllInfo TAR 内的 gzip 和配对 10x MEX，不先展开完整文本。输出：

```text
NEMO_P21_author_20261117/
├── genes/
│   ├── cells.tsv
│   ├── genes.tsv
│   ├── gene_counts.npz
│   ├── isoforms.tsv
│   └── isoform_counts.npz
├── events/
│   ├── cells.tsv
│   ├── genes.tsv
│   ├── groups.json
│   ├── branch_counts.npz
│   ├── branch_offsets.npy
│   └── normalization_factors.tsv
├── conversion_report.json
├── validation.json
└── conversion_status.json
```

配置的 `input.genes`、`input.events` 分别指向以上两个子目录。只有 `conversion_status.json` 为 `complete` 且 `validation.json` 中 `passed=true` 时使用。保留 `conversion_staging.duckdb` 时，重跑可复用已完成的原文件导入事务；原文件未改动是复用前提。中间库不是模型必需输入，本次验证完成后已清理，重新转换会从保留的作者压缩文件重建。

本次已完成结果：**35,317 个细胞、31,053 个基因、3000 HVG、147,447 个局部事件、359,082 个事件分支**，标准输入约 **1.13 GB**。计数、行列轴及模型实际读取检查通过；详细 QC 见转换目录的 `conversion_report.json`、`validation.json`、`model_input_check.json`。对应配置为 `configs/nemo_p21_v0.3.yaml`，本轮没有启动训练。

本次转换采用：

- 细胞键为 `P21_M1_HIPP:barcode` 等，包含样本身份，移除 10x 的 `-1` 后匹配；短读长、作者细胞注释、长读长非冲突 UMI 三者取交集，未匹配细胞数量单列。
- 稳定基因 ID 去掉版本号；全量短读长整数 UMI 用作 G，分母为全基因 library，不能先截成 HVG 后计算 library。
- 按 `sample+barcode+gene+exact_UMI` 去重；同 UMI 有矛盾剪接链时排除并审计，不增加未经验证的编辑距离合并。相同 read ID 有矛盾基因/条码/UMI/链的记录排除。同链 UMI 选实际覆盖跨度最大的读段作为观测，不拼接不存在的分子。
- AllInfo 1-based inclusive 坐标转换为模型 0-based half-open；逐分子核对 intron chain 与 exon gaps。现有数据抽查支持 mm10，不能与主库 mm39 位点直接按数字交集。
- 从观测结构建立 cassette、共享 donor、共享 acceptor 的局部竞争；包含/排除都要求相同两侧各 10 bp 的真实覆盖。默认事件总分子≥10、有效细胞≥3且至少两个分支被观测，不按论文显著性或高变外显子名单挑事件，不强制只有 HVG 可拥有 S。
- 相同 junction chain 的终端覆盖取保守交集，以免 L2 的代表性结构把未覆盖区间当作排除。保留观测链而非宣称全部为全长转录本。
- 3000 HVG 在保留的配对细胞全基因原始计数上按 `seurat_v3` 计算，单独保存；不借用作者聚合结果。模型内部按 CP100000/log1p 归一化，文件保留原始整数。
- 上游 AllInfo 自带作者的比对、条码和结构过滤，不能恢复其已排除分子；本次转换不将这些处理文件称作无过滤原始测序。源流程可含最少结构观测次数等过滤，具体运行门槛未获独立 manifest 确认。
- 当前转换只生成成熟剪接链／事件计数，尚未独立量化 E/I 区域矩阵。`intron chain` 表示已切除缺口，不能直接当作保留内含子或 nascent counts；保留的 AllInfo 比对块可供后续独立区域定量。P21 元数据保留 `trajectory_direction_available=false`，表示来源未提供方向；当前用官方 CytoTRACE 选择每个 DPT 上下文的起点，并记录方向来源。

### Gene / chain 目录

| 文件 | 定义 |
|---|---|
| `cells.tsv` | 唯一 `cell`；可提供 `stage`、`stage_order`、`batch` 和经核验的 Leiden 标签 |
| `genes.tsv` | 唯一稳定 `gene`、`gene_name` |
| `gene_counts.npz` | cell × 全基因整数 UMI；用于完整 library 分母 |
| `isoform_counts.npz` | cell × 实际观测结构的整数分子计数 |
| `isoforms.tsv` | `isoform`、`gene`、JSON `structure`（junction tokens）、JSON `exons`（0-based half-open） |

结构目录中的 isoform 命名不代表所有结构都是全长；模型输出明确为 `observed_splice_chain`。有目标链结构覆盖不等于对所有转录本完整性已验证。

`cells.tsv` 有表头，例：

```tsv
cell	sample_id	barcode	stage	stage_order	batch	donor	region
sampleA:ACGT...	sampleA	ACGT...	day1	0	sampleA	donor1	Hippocampus
```

`cell` 必须唯一，不能只用跨样本重复的 barcode。`stage_order` 是数值先后关系；未知时不随意填出多个阶段顺序。默认 `trajectory.root_method: cytotrace`，也可用 `trajectory.root_cell` 显式指定起点；CytoTRACE 成功选根后可为缺少来源时间方向的上下文运行 DPT。细胞类型注释可另加列，不自动当作训练标签或时间顺序。

`genes/genes.tsv` 有表头，至少 `gene`、`gene_name`，建议增加 `chromosome`、`strand`、`gene_type`；`gene_counts.npz[:,j]` 对应第 j 个数据行。结构表至少包含 `isoform`、`gene`、JSON `structure` 和 JSON `exons`，例如：

```json
{"isoform":"SC_example","gene":"ENSMUSG00000000001","structure":["junction:chr1:100-200:+","junction:chr1:300-500:+"],"exons":[[50,100],[200,300],[500,550]]}
```

实际保存为 TSV 行，两个数组分别编码为单元格里的 JSON 字符串；普通 TSV 引号转义由 pandas/csv 处理。矩阵列严格对应结构表行，负链的 junction token 和 exon 列表按转录方向排列。

### Event 目录

| 文件 | 定义 |
|---|---|
| `cells.tsv` | 与 Gene 输入同一细胞集合；按 ID 对齐，不按位置猜测 |
| `genes.tsv` | 预先确定的 Seurat v3 HVG 稳定 ID 轴；准备阶段不重新挑选 |
| `groups.json` | 每事件包含 `event/gene/gene_name/chromosome/strand/left/right/anchors/branches` |
| `branch_counts.npz` | cell × 连接后的事件分支整数计数 |
| `branch_offsets.npy` | 每个事件分支的起止 offset，长度 events+1 |
| `normalization_factors.tsv` | 有表头的 `cell`、`complete_gene_UMI_library`，建议 `scale_to_100000`；当前 CLI 准备指纹流程要求该文件 |

N=0 是未观测，不当作分支排除。组成模型生成的条件概率只写预测文件，不写回原始 S。不同事件共享分子，目标称为复合计数似然，不假设所有事件独立。

坐标适配支持**小鼠 mm10 与人 hg38**。使用匹配物种的先验数据库，motif 在对应原生 FASTA 重扫，CLIP 使用同一 assembly；检查数据库的 species/primary_assembly 元数据，禁止将小鼠主库当成人源先验。没有隐式 liftover。

`groups.json` 为对象数组，最低结构示例：

```json
[
  {
    "event":"event_001",
    "gene":"ENSMUSG00000000001",
    "gene_name":"Example",
    "chromosome":"chr1",
    "strand":"+",
    "kind":"cassette",
    "left":100,
    "right":500,
    "anchors":[[90,100],[500,510]],
    "branches":[[[100,200],[300,500]],[[100,500]]]
  }
]
```

每个分支是一个有序 junction 边列表，每条边 `[上游外显子结束坐标,下游外显子开始坐标]` 在基因组数字轴上递增；负链的局部分支比较按基因组坐标规范化，链方向另由 `strand` 保存。分支索引没有默认 inclusion 含义，解释时查看结构，不能认定第 0 列必为包含。

若前两个事件分别有 2、3 个分支，则 `branch_offsets.npy=[0,2,5]`，`branch_counts.npz` 有 5 列；第一个事件为 `[:,0:2]`，第二个为 `[:,2:5]`。每细胞该切片的行和为该事件有效分母；没有可判定分子时整段为零，不作为被排除分支。`events/genes.tsv` 仅列 HVG，不要求事件数或事件基因数与 HVG 一一对应。

保存稀疏计数示例：

```python
import numpy as np
from scipy import sparse

sparse.save_npz('gene_counts.npz', sparse.csr_matrix(raw_gene_counts, dtype=np.int64))
sparse.save_npz('branch_counts.npz', sparse.csr_matrix(raw_branch_counts, dtype=np.int64))
np.save('branch_offsets.npy', np.array([0, 2, 5], dtype=np.int64))
```

## 模型流程

1. **G 参考聚类与初始排序**：固定 HVG、CP100000/log1p，PCA→KNN→Leiden；沿用 CausalBridge 的 k=20、使用 20 PCs、resolution=0.5、seed=42。PAGA 默认关闭。各上下文以完整基因计数调用官方 CytoTRACE，取最高分细胞为 DPT 起点；多阶段时在最早阶段内选择。小簇不会删细胞，共享或全局上下文也统一用 **DPT** 排序，不切换 Palantir。不同连通分量和显式分离状态保持边界。
2. **S 局部校准**：以 DPT 为固定初始轴，只修改有可靠剪接证据支持的局部细胞，未支持细胞逐值保持原 DPT 时间。主线为内含子/外显子计数相位证据；可选成熟组成校准。使用带 DPT 偏移惩罚的局部约束、修正幅度限制和阶段边界，不重算第二套全局轨迹。缺少输入时明确记录未校准，不静默改用另一模式。所有事件共享同一份校准时间，全部有支持的基因参与校准，不按目标基因排除特征。
3. **先验**：从同一主库取已有/有标签家族 motif、CLIP、蛋白互作和复合体；按传播轮次预算扩充候选伙伴，并在有 RNA 锚点的伙伴路径上生成间接事件候选。PPI 分数不当作调控方向，co-complex 不冒充直接接触。每目标伙伴预算、阈值、位点来源均保留。
4. **条件响应**：阶段、批次、分支和校准时间背景下，联合拟合少量 RBP 主效应与有互作依据的双因子交互。线性或冻结样条基函数，Dirichlet–multinomial 分支计数目标；浓度为显式共享配置，尚未估计成每事件自由参数。交互对主效应残差化，防止同一变化重复计入。固定状态依赖模式，不要求 Transformer。
5. **虚拟扰动**：在非 log 代理量上乘以 `(1-q)` 再 log1p，冻结背景和其他输入。一次求值联合响应，得到 logit 差分和分支概率差。q=0 与不涉及该事件输入的干预严格不改变输出；联合扰动不把单 KO 简单相加。
6. **结构组成**：只在有实际结构覆盖的细胞/基因上，在合法链 simplex 中重建；公共锚点不足的链不当作排除分支。WT/KO 可行结构优化的相对变化应用到同一 WT 参考，避免 WT 不一致产生虚构 KO 效应。截断结构集合比例、秩和冲突另行输出。
7. **级联**：联合拟合 RBP 总表达与自身局部异构体组成对下游事件的响应。第一轮扰动改变 RBP 自身事件的分支组成后，更新其组成输入，下一轮通过已学习系数传播。各轮同步、相对冻结 WT 更新，直接靶点剂量持续固定，无人工蛋白功能权重或固定时间 lag。

模型拟合使用全部 WT 细胞。输入可缺少单事件计数，此时输出需区分“未观测事件的条件预测”和“已有分子证据”；不称为补全。训练拟合度不代表干预泛化。

### 官方 CytoTRACE 选择 DPT 起点

原版 **CytoTRACE 0.3.1** 已从本地目录复制到 `tools/CytoTRACE_official/`，包含原许可证和未修改的官方 R 源码。`COPY_PROVENANCE.json` 记录原仓库 commit 与逐文件 SHA256。它不使用预训练权重。

```yaml
trajectory:
  root_method: cytotrace
  root_cell: null
  cytotrace:
    source_dir: tools/CytoTRACE_official
    rscript_binary: /home/huangtao/miniconda3/envs/CausalBridge/bin/Rscript
    seed: 42
    enable_fast: true
    ncores: 1
    subsample_size: 1000
```

本机复用已可用的 CausalBridge R 环境，但执行 SpliceFate 下复制的源码，不导入其 Python 包。最小依赖为 R≥4.0、Matrix、HiClimR、ccaPP、nnls；更换机器时调整 `rscript_binary`。Python 包包含 R 传输脚本，官方代码通过配置目录提供，原许可证一并保留。

每个轨迹上下文／连通分量使用**完整原始基因计数**计算评分，不截取为 HVG 或 RBP。取最高有效评分；同分按细胞 ID 固定顺序。有多个已知阶段时，在最早阶段内选择；P21 单阶段直接在该上下文内选择。`root_cell` 在其所在上下文优先。

`enable_fast: true` 使用官方快速计算模式；<3000 个细胞时官方自动关闭，所有细胞均保留。原函数把 `max(mat)<50` 视为 log2 输入，适配器遇到这种低整数计数时先用完整文库转换为 CPM，避免误判；实际输入单位、参数、R 依赖版本与日志记录在输出中。官方评分的 G 平滑不改写原始 G/S 矩阵。

低信息细胞可能没有有限评分：这些细胞保留在模型与 DPT 中，只不参与选根。若整个上下文没有可用评分，则记录不可定向。CytoTRACE 只选起点，细胞顺序仍由 DPT 计算，随后按当前设置做 S 局部校准。

输出 `cytotrace/scores.tsv`、各上下文的 `official.log`、`runtime.txt`，以及 `cells.tsv` 中的评分和 QC。`trajectory_diagnostics.json` 记录 `root_cell`、`root_method`、`root_cytotrace_score`。官方源码、传输脚本和评分文件纳入运行指纹。

仅运行全细胞 CytoTRACE→DPT 参考检查：

```bash
python scripts/check_cytotrace_reference.py --config configs/nemo_p21_v0.3.yaml \
  --output runs/your_new_cytotrace_reference
```

## WT 学习的异构体组成与多轮传播

`prepare` 从候选 RBP 自身的公共锚点事件建立局部异构体分组。每组是同一段真实覆盖下的剪接替代，不把不同终端覆盖直接当成不同完整转录本。默认每个 RBP 保留一个计数有支持的事件、最多四分支；按少数分支总分子数排序，不依据下游响应筛选。支持门槛记录在 `model.isoform_state`，不要求同一 Leiden×样本组各自达到细胞配额。

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

每分支的原始比例为 `p_ik=C_ik/N_i`，丰度输入 `x_i=log1p(CP100000 gene counts)`；组成对比输入为 `x_i × (p_ik-pooled_k) × N_i/(N_i+count_shrinkage)`。K 分支使用 K−1 个对比，最多的分支为参照；因子总表达项仍独立保留。下游模型为这些对比学习有正有负的系数，并可保留与伪时间的交互。不同下游事件有不同的权重，不将其解释为该异构体全局固定的蛋白功能分数。

仅在上游组成和下游事件都实际观测到的细胞上检查可拟合性、确定组成基函数。原有全部 WT 细胞仍参与事件拟合：上游未观测时，该细胞的组成项被观测门控设为零，另含观测门控背景项，不给它补一个比例。一个因子的自身事件不使用该事件定义的自身组成作为预测输入；这个规则只作用于响应模型，不排除其参与共享伪时间校准。无足够组成信息时保留总表达模型。

`fit` 导出 `isoform_response_weights.tsv`，每项包含上游 RBP、上游事件／分支、下游事件／分支、标准化尺度、联合有效细胞数和拟合系数。`fit_summary.json` 记录真正使用组成输入的事件数量。

模拟先在非 log 表达量上执行 KO/KD，再同步调整该因子的表达加权组成项。每轮使用冻结模型计算相对 WT 的事件 log 概率变化，将该变化作用到上游有观测的组成；总表达除直接剂量外保持不变。组成更新后再次计算下游响应。不会逐轮叠加同一效应，不执行人工功能权重，也不需要 domain/CDS 注释。无新组成输入变化、没有支持的组成项或达到轮数预算时停止。

`simulation.json` 保留每轮改变的表达输入、组成输入、上游事件及更新细胞数；`isoform_state_diagnostics` 明确区分“有组成输入”“源事件已拟合”和“被下游模型实际使用”。`--bootstrap` 对每个 WT 重采样拟合全部参与事件并重新运行组成级联，输出稳定性区间。所有事件共享已冻结的一维校准时间。

旧的 `input.functional_feedback` 人工权重接口已停用；传入非空旧路径会明确要求迁移，不再读取它。

## scVelo 与长读长动力学

本模块的主要目标是借鉴**计数关系→方向证据→细胞时间校准**的逻辑，自行实现所需部分；完整 scVelo API 适配仅保留作可选比较。

### 独立内含子计数校准

`src/splicefate/intron_time.py` 无 scVelo/Regionvelocity 依赖、不复制其源码。`input.intron_counts` 指向以下目录：

| 文件 | 定义 |
|---|---|
| `exon_counts.npz`、`intron_counts.npz` | 同一 cell×gene 轴，整数 UMI 分子或明确的 molecule-region incidence |
| `cells.tsv` | `cell`、`normalization_denominator`（全量库的共享分母） |
| `genes.tsv` | 稳定 `gene` ID，不重复 |
| `quantification.json` | assembly、annotation_source、classification_validated、count_unit、intron_definition、exon_definition、ambiguity_policy、normalization_denominator_definition |

每个基因在有支持的高计数细胞中估计 `I ≈ slope · E`，以 `I - slope · E` 提供相位偏离方向；与候选邻近细胞的外显子状态差比较，计算覆盖加权方向一致性。随后只对这些局部细胞求小幅 DPT 修正：约束相邻次序与方向的一致性，同时惩罚偏离原 DPT，固定根及阶段边界。没有证据的细胞保持原值，不通过重新计算全局最短路带动其他细胞。该方向是辅助证据，不是完整动力学速率，也不是内含子总数年龄排序。

不平滑或补写原细胞计数，不以 retained intron 的存在一概判定 nascent RNA。仅对两类区域都有分子支持的基因进行局部计数关系判断；保留区域定义、歧义处理和定量来源。长读长每分子包含的内含子数可以贡献 incidence 统计，但不能只凭单一内含子总数提供方向，仍需对应外显子量与背景。

**当前 NEMO P21 数据的 E/I 情况**：配对 MEX 是每基因总 UMI，AllInfo 的 intron chain 是剪接后基因组缺口；两者都不能直接作为独立的未剪接内含子计数。2026-11-20 可行性抽查四份 AllInfo 各前 10,000 行，发现 830、1,011、807、921 行的比对块覆盖了所属基因“所有参考转录本外显子并集以外”的内部区域≥10bp。这是候选 intronic coverage，说明不能断言原文件完全没有内含子信号；抽样不是随机覆盖率估计。正式区域矩阵还需排除其他基因／外显子的歧义、逐分子 UMI 去重和明确区域定义，且不把区域保留自动标成 nascent。该独立定量尚未执行；`intron_counts: null` 时继续使用 CytoTRACE 定向的 G-only DPT，不生成假的 E/I 值。已有 S 事件计数照常用于响应拟合和组成级联。

`trajectory.intron_phase` 显式控制支持门槛、参考分位数和辅助权重。已知阶段只提供宏观锚点，分支内细化，所有事件使用同一份冻结时间。该初版使用相位残差，而非完整 dynamical EM；未来可以独立比较完整动力学，不将二者混报。

### 可选完整动力学比较接口

存在专门的长读长方法 **Region Velocity**，见 [NAR 2023 论文](https://doi.org/10.1093/nar/gkad969) 与 [源码](https://github.com/Dekayzc/Regionvelocity)。其输入为外显子/内含子区域计数，并可使用 RNA 分子数计算 theta，不是成熟异构体包含/跳过计数。

可执行的可选 scVelo 适配器位于 `splicefate.velocity`：

```bash
PYTHONPATH=src /path/to/scvelo/python -m splicefate.velocity \
  --h5ad validated_precursor_mature_counts.h5ad \
  --prepared runs/your_WT_reference \
  --output velocity_calibration.tsv
```

输入 `layers['spliced']`、`layers['unspliced']` 必须是真正基因级整数 UMI，并在 `uns['velocity_input_qc']` 中记录 `classification_validated=true`、`unit=umi`、assembly、annotation_source、quantification_method、precursor_definition、mature_definition。IR、成熟 AS 分支、无剪接 junction 的单外显子 reads 不能直接替代该分类。

适配器调用 scVelo `recover_dynamics(use_raw=True, fit_connected_states=False)`、dynamical velocity、velocity graph 和 latent time。使用动力学 likelihood≥0.1，不额外套用稳态 R² 门槛；不调用 `pp.moments` 替换观测。保留的图/latent-time 正则与原始分子观测分开。每个轨迹上下文使用全部合格基因拟合一次，输出共享的 `latent_time` 与 `qc_pass` 列；失败上下文保留空值及原因，全部失败时返回非零退出码。

可将输出路径填入 `input.velocity_times`，以 `trajectory.velocity_weight` 辅助状态方向；qc 未通过的细胞退回原 G+S。scVelo 是可选依赖，不安装于主 package 的必需依赖链。Region Velocity 源码已审查，尚未提供 Region Velocity 原算法运行包装器或重写其 GPL-3 实现。

## 输出与审查

| 文件 | 内容 |
|---|---|
| `provenance.json` | 输入、代码、依赖、配置和准备产物指纹 |
| `prepared_arrays.npz` | `times` 为长度等于细胞数的一维共享伪时间；拟合、bootstrap 与所有事件模拟统一读取 |
| `cells.tsv`、`trajectory_diagnostics.json` | 参考聚类、上下文、G/GS 时间、方向和批次提示 |
| `priors.json` | 位点和伙伴候选、证据来源、有效因子轴 |
| `fit_diagnostics.tsv`、`all_event_status.tsv` | 拟合状态和所有事件未执行原因 |
| `isoform_response_weights.tsv` | WT 学习的上游局部异构体组成对下游事件的有符号响应系数 |
| `models.json`、`structures.json` | 冻结模型、输入支持、结构及 eligibility |
| `simulations/*/events.tsv.gz` | 逐细胞事件分支概率、差分、分母和外推/识别标记 |
| `simulations/*/chains.tsv.gz` | 有观测支持的链组成变化及保留比例 |
| `event_stage_summary.tsv` | 同一细胞集合汇总，携带靶点与剂量 |
| `simulation.json` | 轮数、停止原因、解释范围和不确定性状态 |

未运行 bootstrap 时不填造区间。bootstrap 是 WT 细胞按 stage×batch 重采样稳定性，不等于跨供体误差或真实 KO 因果置信区间。外部评价表必须同时匹配 `targets,strength,event,cell_stage,branch_index`，避免错误比较不同剂量或靶点。

## 参考源码与许可

- [CausalBridge](https://github.com/Senki-zip/CausalBridge_en)，固定 `821c0314850fe12fde9d017d41c577bc4c7efd59`；借鉴先验路径、冻结 WT 差分、小型响应与来源审计逻辑。本实现不是直接调用其 ATAC 模型。
- [scVelo](https://github.com/theislab/scvelo)，审查源码 `f63c0e70596ced2f1bee8cf07e8ab66037cf86b2`；适配器 API 检查版本 0.3.4。
- [Regionvelocity](https://github.com/Dekayzc/Regionvelocity)，审查源码 `533483156993300e47e9aeb63662977388771cec`，版本 1.3、GPL-3。只读参考，不将 GPL 算法源码拷入 Python 包。

第三方源码克隆在 `tools/`，保留原许可。本轮没有从这些工具复制数值算法代码，scVelo 通过其公开 API 可选调用。
