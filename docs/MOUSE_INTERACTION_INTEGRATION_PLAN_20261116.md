# 小鼠剪接因子互作数据：下载、覆盖审查与整合方案

日期：2026-11-16。新增数据目录：`data/external/rbp_interactions_20261116/`。

实施更新：用户确认后，关系数据已正式导入小鼠主库。当前使用入口和实际表结构见 `data/external/mouse_rbp_unified_20261116/INTERACTIONS.md`，验收见同目录 `interaction_integration.json`、`validation.json`。下文保留下载预审时的计划与口径；人类阶段投射、逐篇 WT/RNase 审查仍待后续工作。

下载已完成：11 个核心导出文件，共 320,057,621 字节；附加保存 SpliceosomeDB 页面与人类复合体阶段比较矩阵。ZIP CRC、gzip 完整性、文件 SHA256、表结构及初步关联检查通过。

本次建立了独立的 `audit/interaction_preview.duckdb` 预整合库和逐因子覆盖表。原 `mouse_rbp_unified_20261116/mouse_rbp.duckdb` 以只读方式访问，尚未修改其表或接入训练。该预整合库用于评估关系层扩充，不把网络连接当作 motif 或因果调控边。

## 1. 实际取得的数据

| 来源 | 本地数据 | 规模与解释 |
|---|---|---|
| IntAct 当前版 | `raw/intact/mouse.zip` | 206,564 条记录；其中双方蛋白来源物种均为小鼠的记录 34,470 条。筛选正向蛋白关联/结合类型后为 30,684 条证据、21,107 个无向源蛋白 ID 对。多证据记录不计为独立蛋白对。 |
| IntAct 阴性结果 | `raw/intact/mouse_negative.txt` | 单独保存，不并入正向互作网络。 |
| STRING 12.5 | 小鼠 physical full、functional detailed、蛋白信息、别名 | 物理网络去掉 AB/BA 后 390,379 个蛋白对；功能网络 6,753,475 个蛋白对。保留所有通道和转移评分，不能把功能关联全部归为物理接触。 |
| Complex Portal | 小鼠 `10090.tsv`、人类 `9606.tsv` | 小鼠 753 个复合体、人类 2,533 个。小鼠条目中仅 78 个明确标为物理互作证据，其余分别属于背景知识、orthology、homology、paralogy 或混合物种实验，标签全部保留。 |
| SpliceosomeDB | 小鼠/人类蛋白 CSV、全物种实验 CSV、人类阶段矩阵 HTML | 小鼠 340 条蛋白目录；人类 1,110 条，137 个条目有阶段归属，共 430 条归属；实验目录 168 条。小鼠目录的阶段字段为空。 |
| 已有 CORUM | 从当前统一库读取原始支持记录 | 复用 671 条 Mouse 复合体注释，避免重复下载和重复计算来源。 |

物种字段记录分子来源；实验宿主细胞、组织与条件单独记录。宿主、方法、实验角色、扩展方式等 MITAB 字段完整保存在 `record_json`，正式导入时需要结构化。

小鼠 Complex Portal 证据分布：物理 78；背景知识推断 273；orthology 318；homology 62；paralogy 8；混合物种实验 14。各条目按物理实验、推断等来源类型分别统计。

## 2. 当前因子覆盖

沿用去别名后的 406 个 GO 宽泛剪接候选，其中 104 个有原始或家族扩展 motif，302 个没有当前库的 motif 关联。“没有收录”不等于该蛋白没有序列偏好。

| 关系证据口径 | 406 个候选中覆盖 | 302 个无 motif 候选中覆盖 |
|---|---:|---:|
| IntAct 小鼠–小鼠正向蛋白关联/结合类型 | 260 | 190 |
| 其中明确标为 direct interaction | 10 | 9 |
| STRING physical 中非转移实验通道有分数 | 333 | 245 |
| 同一通道分数至少 400（敏感性统计） | 321 | 237 |
| 同一通道分数至少 700（敏感性统计） | 69 | 45 |
| Complex Portal 物理证据或 CORUM Mouse 复合体 | 24 | 15 |
| IntAct、STRING 非转移实验通道及上述复合体的并集 | 339 | **249** |
| STRING 非转移实验通道直接连到一个有 motif 的因子 | 175 | **105** |

400/700 用于评估覆盖对评分口径的敏感性。STRING 分数按提供方的互作证据评分保存，组织与实验背景单独记录；跨平台来源重叠保留追溯信息。

因此，302 个无 motif 候选中，249 个（82.5%）至少能补充一类上述关系先验。它们与 104 个有 motif 因子合计，353/406（86.9%）有 motif 或此口径关系信息。该比例按因子级先验覆盖统计。

Ddx5、Hnrnpm、Prpf8、Sf3b1、U2af1 等都能取得此类互作信息，尽管当前 motif 表没有相应关联。逐因子明细位于 `audit/splicing_factor_interaction_coverage.tsv`。

## 3. 建议的正式整合结构

保留现有因子、motif、位点和 CLIP 层，在同一数据库中增加独立的关系表。

| 拟新增表 | 关键字段 | 用途 |
|---|---|---|
| `interaction_sources` | 来源、版本、URL、快照日期、文件校验值 | 可复现的输入登记 |
| `interaction_identity_map` | 来源 ID、稳定基因 ID、映射方法、歧义集合、原物种 | 复用规范化小鼠身份；歧义不强制归并 |
| `protein_interaction_evidence` | 两端因子、类型、检测方法、互作 accession、PMID、分子物种、宿主、细胞/组织、RNA 依赖状态、阴性标记 | 保存每条实验或来源证据，不丢失检测背景 |
| `protein_interaction_pairs` | 规范化无向蛋白对、支持证据引用 | 去除 AB/BA、来源重复视图；不直接累加数据库分数 |
| `complexes` | 来源复合体 ID、名称、物种、证据类型、参考文献 | 保留完整复合体及推断状态 |
| `complex_memberships` | 复合体、成员、化学计量、来源证据 | 多成员关系；不自动扩成全部两两直接接触 |
| `spliceosome_stage_memberships` | 因子、复合体/阶段、原物种、映射证据、文献 | 保存阶段归属；人类与小鼠分开 |
| `interaction_prior_views` | 上述关系的可筛选集合 | 按证据类型选取先验，支持回溯与消融 |

现有预整合库已提供这些表的大部分原型，但不替代正式身份、实验来源和证据去重审查。

### 证据标签

至少区分：`direct_physical`、`physical_association`、`association_unspecified`、`co_complex`、`proximity`、`colocalization`、`functional_association`、`orthology_inferred`、`homology_inferred`、`paralogy_inferred`、`background_knowledge_inferred`、`mixed_species_experiment`。

RNA 依赖性使用 `RNA_dependent / RNase_resistant / unknown`，只有存在明确实验才填前两类；未报告 RNase 处理不等于 RNA 非依赖。检测到 RNase 抗性也不单独证明两个蛋白直接接触。作用方向和效应大小默认未知，不根据无向 PPI 填正负号。

### 去重与身份

1. 复用当前物种内稳定基因 ID 和已记录的别名规范化。STRING 映射优先稳定基因 ID，再用唯一 UniProt，最后用唯一精确符号。无法解析的来源记录保留。
2. 同一无向蛋白对存一份，实验支持另存。优先按原始 interaction accession/实验 accession、方法与背景识别同一证据；PMID 相同用于提示复用，不能把一篇文章中的不同实验全部合为一个。
3. IntAct、Complex Portal、CORUM、STRING 可能复用同一实验，不进行独立性累加。对无法确认是否相同的证据保留来源标签。
4. 人类复合体和阶段先保留原物种。后续若引入可靠一对一 ortholog，另存 `orthology_inferred` 关联，不覆盖原始小鼠证据；RBP2GO 的广义 homolog 列表不能自动当成一对一 ortholog。

## 4. 怎样扩充模型可用信息

建议先将新库用于识别“因子—复合体/模块—剪接过程”关系，并保留两条入口：

- 有 motif/CLIP 的因子：保留当前因子–RNA 位点候选。
- 无 motif 的因子：由有来源的互作或复合体归属提供模块约束，之后结合本数据中剪接事件的变化学习关联。

105 个无 motif 因子能在 STRING 非转移实验通道中连接到有 motif 的伙伴，可以优先检查，但这只是**候选机制路径**。不能把伙伴的全部 motif 复制给它，也不能宣称它调控伙伴的所有靶事件。核心剪接组分可能影响识别、组装、催化或重塑，其效应应允许作用于模块和多个事件。

这一关系层既可作为稀疏回归的结构约束，也可用于模块评分或图模型；是否采用图神经网络属于后续模型比较问题，不由数据库形式决定。

## 5. 建议的实施顺序与验收

1. **先正式导入小鼠身份明确的实验关系和复合体成员**，保留其他类型供筛选，提供按因子查询互作和复合体的接口。
2. **再做来源去重和人工机制抽查**：先检查 Ddx5、Hnrnpm、Sf3b1、Prpf8；核对实验方法、物种、共复合体与直接接触的区别。
3. **随后评估人类阶段信息的有标签投射**，依赖可靠 ortholog 映射。不能用人类字段填补小鼠而不留标签。
4. 建立最终证据覆盖表和 provenance 查询，确认未知 RNA 依赖性、阴性结果、歧义 ID、同源推断和被复用实验可单独筛选。

当前不把扰动响应数据加入 WT 训练先验，不改训练参数，不以这些覆盖数字替代独立 KO/KD 验证。

## 6. 本轮未取得或未纳入的资料

- Hegele 2012 人类剪接体原始网络：公开 IntAct 文献查询返回空结果，不能据此说已下载论文全部网络。仍可进一步寻找作者原始补充表。
- RINE 为酵母资源，本轮小鼠证据层未导入，也未投射为小鼠网络。
- Papasaikas 等功能网络来自扰动实验，本轮未作为 WT 先验导入。
- 小鼠剪接阶段的原生证据在下载的 SpliceosomeDB 条目中为空，当前覆盖为零。

完整下载状态见 `download_status.json`；逐文件来源见 `download_manifest.json`；预审结果见 `audit/coverage_summary.json`。
