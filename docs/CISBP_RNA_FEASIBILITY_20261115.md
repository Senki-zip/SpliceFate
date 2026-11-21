# CISBP-RNA 对 SpliceFate 的适用性核查（2026-11-15）

CISBP 提供模型所需的 motif 候选先验；WT 数据用于拟合条件响应，独立扰动数据用于评价。

## 数据与统计口径

从官方网站按物种下载 build 2.00 的 RBP Information 与 PFMs；保留 ZIP、下载 URL、SHA256 和逐矩阵 QC。统计按唯一 RBP_ID，D/I/N 使用数据库自身标签，分别表示直接测量、推断、无 motif；覆盖统计以下载包中的唯一 RBP_ID 为单位。

| 项目 | 小鼠 | 人 |
|---|---:|---:|
| RBP 条目 | 424 | 465 |
| 有 motif 记录 | 209 | 271 |
| 默认关联中有非空可解析矩阵的 RBP | 200 | 265 |
| 全部关联中有非空可解析矩阵的 RBP | 208 | 267 |
| 下载包中的空矩阵文件 | 38 | 42 |
| 元数据引用但包内缺失的 motif ID | 5 | 132 |

小鼠状态为 D=7、I=202、N=215；人状态为 D=189、I=82、N=194。默认与替代 motif 分别记录，空或缺失文件按下载包状态统计。小鼠 Nova2 在该包的关联矩阵均不可用，可另从原论文或其他 motif 资源补充。

## 与现有小鼠数据的交集

使用完整的 161×32285 基因计数矩阵，按 Ensembl gene ID 精确映射。209 个有 motif 记录的 RBP 中，199 个进入当前基因注释，192 个至少一个细胞有观测计数，174 个至少十个细胞有计数，25 个位于旧版 3000 HVG。统计范围为全体 RBP；十个细胞是覆盖统计阈值。具有可用关联矩阵且至少一个细胞检出的有 191 个。

| 因子 | 数据库证据 | 检出细胞/161 | 旧版 HVG |
|---|---|---:|---|
| Ptbp1 | I | 149 | 否 |
| Ptbp2 | I | 82 | 是 |
| Rbfox2 | I | 95 | 否 |
| Celf1 | I | 138 | 否 |
| Esrp1 | I | 95 | 是 |
| Srsf1 | I | 154 | 否 |
| Srsf2 | I | 159 | 否 |
| Srsf3 | I | 161 | 否 |
| Srsf4 | I | 103 | 否 |
| Srsf5 | I | 157 | 否 |
| Srsf6 | I | 156 | 否 |
| Tra2a | I | 143 | 否 |
| Tra2b | I | 158 | 否 |
| Nova2 | D | 1 | 否 |
| U2af1 | N | 159 | 否 |

表达覆盖按观测 UMI 统计。RBP 输入从完整基因矩阵中单独提取。

## 能力与缺口

- 可提供：因子 ID、物种、结合结构域、结合偏好矩阵、实验/推断来源、部分推断方法说明。配合参考基因组及事件坐标，可定位候选顺式元件。
- 尚不提供：目标细胞中实测占据位置、每条 RBP–event 边的促进/抑制方向、KO 剂量响应、蛋白有效活性、目标系统真实扰动效应或命运终点。
- 相同 motif 不能唯一归因：小鼠 Srsf4/Srsf6 默认共享 M00815_2.00；Tra2a/Tra2b 共享 M00826_2.00。需要结合位置、非扰动 CLIP、状态条件和独立验证，由响应模型估计作用方向。
- 候选位置应在转录方向正确的 pre-mRNA 上扫描，包含事件相关外显子及内含子窗口，命中记录作为候选先验。
- WT-only 制度可使用非扰动结合先验。由 KD 学得的 KATMAP 等有方向权重须归入独立迁移模式；盲测效应不得用于选边或调参。

## 开发判断

可以以 CISBP 为主要 motif 来源启动小规模先验和状态条件模型开发。接下来应对候选因子做表达变异/共线性检查、与事件序列和物种匹配 CLIP 求交，并锁定独立扰动验证。候选表达较充分的因子包括 Ptbp1、Rbfox2、Celf1、Esrp1、Srsf1/2/3；这里只按 motif 与表达给出候选，尚未按 CLIP、变异度或验证可得性完成排序。

本次核查范围为 motif 文件、因子身份与表达覆盖。

## 来源和可复核文件

- [CISBP-RNA](https://cisbp-rna.ccbr.utoronto.ca/)
- [官方帮助与 Direct/Inferred 定义](https://cisbp-rna.ccbr.utoronto.ca/help.html)
- [官方下载](https://cisbp-rna.ccbr.utoronto.ca/bulk.php)
- [2025 年资源论文](https://doi.org/10.1038/s41587-025-02733-6)
- 数据目录：`data/external/cisbp_rna_2.00/`；含 `download_manifest.json`、`audit_summary.json`、`matrix_qc.tsv`、`mouse_factor_coverage.tsv`、`mouse_selected_factors.tsv`、`shared_default_motifs.tsv` 和人/小鼠原始 ZIP。
