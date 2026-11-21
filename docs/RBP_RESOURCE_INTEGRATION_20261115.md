# RBP 数据下载与整合（2026-11-15）

已整合 CISBP-RNA、ATtRACT、RBP2GO、ENCORI 和 POSTAR3 的公开可取得数据。
POSTAR3 的整库位点坐标仍需向提供方申请，本次未提交申请；已完成公开逐因子查询。

| 来源 | 本次取得的数据 |
|---|---|
| CISBP-RNA 2.00 | 全物种 RBP/全部关联 motif、29,442 份可解析矩阵、SQL 备份 |
| ATtRACT | 全量 4,882 条结合序列注释和矩阵 |
| RBP2GO 2.0 | 13 个物种的 29,215 条 RBP 和 147,681 条背景蛋白记录；同源、GO、结构域等公开导出 |
| ENCORI | 2,725 个 CLIP 实验的 BED 位点，12 个 RBP–RNA 注释包，公开 HOMER motif 候选 |
| POSTAR3 | 七个物种 1,434 次公开查询，351 个因子有非空结合结果；整库坐标待申请 |

统一后的主要规模：276,233 个精确矩阵组、406,866 条来源可追溯的因子–motif 连接、6,119,415 条因子–基因/类型结合注释。这些连接记录结合证据和来源。合并因子表包含背景蛋白及未解析标识，不能把总行数解释为 RBP 数量。

小鼠中按 GO 描述筛得 408 个宽泛的剪接相关候选，其中 393 个与当前基因注释精确匹配。此列表需结合表达、细胞系统和扰动验证进一步筛选，未修改模型或基因选择参数。

## 整合与验证

- 因子按物种和基因 ID 匹配，保留别名与映射方法；歧义项不强行合并，跨物种不合并。
- 矩阵保留实验/推断证据与原始注释；空、缺失矩阵单独记录。同一数值矩阵去重，来源连接保留。
- 位点保留原始 genome build、链方向和 BED 的 0-based half-open 坐标；不做跨物种坐标转换。
- ENCORI 记录原始实验 accession、相同文件内容和重复区间；POSTAR 聚合结果缺少实验 ID，跨来源重叠按支持记录保存。
- ENCORI 原始包中的 KD/CRISPR 表达响应列保留在原始归档内，未导入 WT 结合先验表。
- 文件 SHA256、ZIP CRC、矩阵、BED 行、SQLite 完整性和关联标识均已核查。验证通过：True。
- CISBP 全物种导出补回了早期物种包中缺少的部分矩阵。

## 使用入口

主数据库：`data/external/rbp_atlas_20261115/integrated/rbp_atlas.sqlite`

```bash
python scripts/query_rbp_atlas.py --gene Ptbp1 --species Mus_musculus
```

数据库表和使用边界见 `data/external/rbp_atlas_20261115/README.md`。
来源状态见 `source_status.tsv`，逐文件清单见 `download_manifest.tsv`，验证结果见 `validation.json`。
小鼠候选覆盖见 `integrated/mouse_splicing_factor_coverage.tsv`。

本次下载整合只涉及 RBP/剪接相关资源，不包含各平台无关的 miRNA、癌症生存分析或原始测序 reads 整库镜像。
