# GSE295932 数据选择与整数 UMI 簇转换

更新日期：2026-11-21。当前开发决策：G、剪接链和局部事件分支统一使用作者发布表的 **`cluster` 整数列**，不使用 `count` 加权列，不使用 `size` reads 列。配置为 `model.count_unit: integer_umi`。本次仅重新转换及检查，不启动建模。

## 1. 细胞和文库选择

用于拟合的是 CRISPR 实验的全转录组非靶向对照，按作者 `barcode → target` 分组表中 `target=Nontarget` 筛选：

| 条件 | 样本 | 有计数细胞 | 原文件记录行数 |
|---|---|---:|---:|
| 未刺激 | GSM8961552 | 387 | 1,806,477 |
| 刺激后 | GSM8961554 | 662 | 3,078,071 |

两组共 1,049 个细胞。原刺激分组表有 663 个 Nontarget 条码，其中 1 个在全转录组计数文件中没有记录，不能补成有观测细胞。细胞 ID 使用 `GSM:barcode`，不同文库同 barcode 不合并。Nontarget 保留原标签，表示非靶向 CRISPR 对照，不是逐细胞验证的编辑状态。

不纳入 KO 细胞；不纳入 GSM8961553/55 的靶向富集文库；不纳入 GSM8961556/57 的单独 WT ONT/PacBio 基准文库。后两者来自同一份 10x cDNA 的平台分样，去重后有 5,726 个细胞，但基因覆盖由靶向富集决定。对应分离及合并记录保留在 `test/GSE295932/WT/`。

## 2. 三个数值字段的含义与选择依据

核对论文方法、发布的 GEO 表、[作者 LongcellPre v1.0.0 代码](https://github.com/yuntianf/LongcellPre/blob/v1.0.0/R/reads_filter.R)及 Zenodo Fig. 6 对象后，字段对应如下：

| 字段 | 作者代码中的汇总 | 统计单位 | 本次使用 |
|---|---|---|---|
| `size` | `sum(size)` | 支持 reads，包含 PCR 扩增 | 不使用 |
| `cluster` | `n()` | 分配到校正结构的 UMI 簇数量，整数；小簇散射过滤之前 | **G、链、事件均使用** |
| `count` | `sum(weight)` | UMI 小簇／散射校正后的加权计数，发布表存在分数 | 不使用，保留对照审计 |

`cluster` 不是细胞聚类编号。它由作者逐 UMI 簇校正结构后按细胞／isoform 汇总；原始 UMI 字符串没有保存在这张汇总表中。它可提供整数、分子簇层面的输入，但不是未经处理的 reads，也不是已独立验证的真实分子 whitelist。

| 条件 | 使用的 `cluster` 总数 | 旧 `count` 总权重 | 差额 |
|---|---:|---:|---:|
| 未刺激 | 2,087,100 | 1,986,503 | 100,597 |
| 刺激后 | 3,599,481 | 3,439,708 | 159,773 |
| 合计 | **5,686,581** | 5,426,211 | **260,370** |

差额来自两种计数定义，不能解释成新增测序或恢复了所有真实分子。特别是 `count=0, cluster>0` 的发布行，在本次整数版本中按其 cluster 数保留，不沿用加权版的正 count 筛选。逐样本相关行数、source SHA256、选择字段及对照权重和均写入 `conversion_report.json`。

## 3. 为什么不用加权 count，也不直接四舍五入

本项目原计数目标是整数分子／分子簇的分支计数。之前的转换保留了 `count` 的分数，并增加 gamma 扩展 DM 伪似然接口；这可以处理估计权重，但改变了计数假设。本次恢复整数输入，通过选择作者已发布的整数 cluster，而不是对 count 取整或按 reads 复制分子。

旧加权版保留在 `data/processed/GSE295932_Nontarget_hg38/`，仅作历史复核；当前配置改指向 `data/processed/GSE295932_Nontarget_hg38_cluster_20261121/`。没有覆盖旧矩阵，两个版本的计数定义不可混用。早期含整数截断错误的尝试保留在带 `superseded_integer_helper` 后缀的目录，验收已标记无效，也不作为当前输入。

本次重新统计结果：

| 指标 | 历史加权 count 版本 | 当前整数 cluster 版本 |
|---|---:|---:|
| 细胞 | 1,049 | 1,049 |
| 基因 | 16,699 | 16,699 |
| 观测剪接链 | 217,053 | 217,053 |
| 局部事件 | 28,261 | **28,266** |
| 事件分支 | 60,030 | **60,040** |
| G 总计数 | 5,426,211 | **5,686,581** |
| 观测剪接链总计数 | 约 4,029,612.89 | **4,237,194** |

细胞及基因 ID 轴完全一致；两版各 3,000 HVG 中重叠 2,878 个，122 个被替换。改变来自完整 G library 的计数口径，而非手动改选基因方法。当前事件分支计数总和 13,289,750，可重复使用同一个来源簇，不当作独立总分子数。全部矩阵均为 int64、分数项为 0。新增事件是相同发现规则下达到计数保留门槛的结果，没有改变 anchor 或事件过滤参数。

`cluster` 在散射／小簇过滤前，可能包含错误分散产生的小 UMI 簇。其整数性质不代表已解决这些噪声；后续结果可与加权版做敏感性比较。当前步骤不自行重做 UMI 聚类，不假定汇总后的整数就能还原逐 UMI 的归属。

## 4. G、结构链和事件的生成

转换脚本：`scripts/convert_gse295932.py`。结构标识、事件发现、公共锚点计数复用 `splicefate.convert`。

1. 校验 cluster 有限、非负且为整数；小数、负数、缺失值直接拒绝，禁止取整补救。
2. G 按 `sample + cell + gene` 对所有作者分配到该基因的 cluster 求和。没有配对短读长矩阵，因此这是**长读长基因总 UMI 簇数**。归一化使用完整 G library，CP100000/log1p；HVG 从本版本 G 重新计算 Seurat v3 3,000 个，不复用旧版列表。
3. 参考为本地 `refdata-cellranger-arc-GRCh38-2024-A`、GENCODE v44。作者外显子区间从 1-based inclusive 转为 0-based half-open，gene ID 去版本号。
4. 同一 gene＋junction chain 合并整数 cluster；终端覆盖取同链记录的保守交集。单块结构参与 G，但不伪装为多外显子剪接替代；不声明完整转录本。
5. 发现观测到的 cassette、共享 donor／acceptor 局部竞争，每个分子簇对应的作者结构必须同时覆盖两侧公共 10 bp 锚点。事件保留总 cluster≥10、有效细胞≥3、至少两分支有观测，不使用论文显著性名单或 KO 响应筛选。
6. G、链和分支矩阵均保存 int64。不同事件可复用同一个来源分子簇，所以事件计数总和不应与 G 总和直接相等；逐事件核对来源贡献和分支总量。
7. 参考未解析的基因继续留在 G 并单列审计；缺少可解释坐标的基因不生成 S。S 缺失不当作分支排除或填补。

## 5. 模型及状态设置

当前人源配置为 `configs/gse295932_nontarget_hg38.yaml`，默认启动脚本读取该配置：

- `clustering.enabled: false`，不运行 Leiden；PCA、邻接图、CytoTRACE 选根和 DPT 保留。
- `priors.targets: [ALL]`；实际模拟靶点 RBFOX2，强度 0、0.5、1。
- 刺激／未刺激保留各自上下文，不把刺激设为必然发育时间顺序。
- `model.count_unit: integer_umi`，使用整数分支 DM 复合计数目标；兼容旧加权模式的代码仅供历史数据选择，不由当前配置启用。
- `intron_counts: null`：cluster 提供的是剪接结构计数，不是未切除内含子区域计数；独立 E/I 相位矩阵仍未量化。

## 6. 验收与复现

```bash
conda activate SpliceFate
python scripts/convert_gse295932.py
python scripts/run_splicefate.py --check
python -m pytest -q
```

默认输出为新的带 cluster 日期目录，完成目录拒绝覆盖。`quantification.json` 固定计数来源；`conversion_report.json` 记录逐样本选择字段和来源；`validation.json` 验证所有矩阵整数性、全部细胞、完整库分母、G 总 cluster 守恒、链统计守恒、各事件总量与唯一标识；另用逐结构标量统计抽查分支计数，重新读取 SHA256 和真实模型输入接口。

论文：[Single cell and spatial alternative splicing analysis with Nanopore long read sequencing](https://doi.org/10.1038/s41467-025-60902-2)。图源数据：[Zenodo 15320816](https://doi.org/10.5281/zenodo.15320816)。本地文献核查与整数列审计：`test/GSE295932/metadata/integer_UMI_literature_findings.json`、`integer_UMI_cluster_field_audit.json`。
