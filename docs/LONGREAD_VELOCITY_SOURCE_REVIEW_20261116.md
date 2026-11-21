# 长读长 velocity 与剪接辅助伪时间源码核查

本轮在实现途中按用户要求检索和阅读源码，随后按“先不管当前数据集”的要求，只继续独立合成验证，不以已有数据的样本量或输入缺口否定通用模型路线。

用户进一步明确：借鉴其利用内含子数量关系校准时间的实现逻辑，不需要全套源码调用。当前主线已改为独立 `splicefate.intron_time`：基因级 I/E 参考关系、相位偏离、候选有向图和冻结时间校准。它采用辅助相位方向，不等同于复现 scVelo EM 或 Region Velocity 全算法；独立分子/区域计数、归一化分母与共享校准时间均可检查。scVelo 完整接口保留作比较，不作为主线依赖。

## 有直接适配长读长的方法

**Region Velocity**，论文 *Improving the RNA velocity approach with single-cell RNA lifecycle (nascent, mature and degrading RNAs) sequencing technologies*，NAR 2023：[论文](https://doi.org/10.1093/nar/gkad969)，[GitHub](https://github.com/Dekayzc/Regionvelocity)。支持 ONT/FLOUR-seq 和 PacBio/HIT-scISOseq V2，属于直接相关的参考实现。

固定源码 `533483156993300e47e9aeb63662977388771cec`，本地 `tools/Regionvelocity_reference`。已读 README、DESCRIPTION 及 `R/momentum_routines.R`：

- `gene.region.velocity.estimates`（约 708 行）：输入 exon/intron 区域矩阵，共有基因和细胞轴对齐；可用 RNA 分子矩阵计算 theta。
- theta 分支（约 894 行）：外显子、内含子计数分别除以 RNA counts，体现 region 数与 molecule 数的区别。
- `gene.EM.velocity.estimates`（约 1032 行）：从初始区域响应迭代求解并输出投射状态。
- `single.gene.EM.inference`（约 2899 行）：逐基因诱导/衰减时间和参数迭代。

该包默认 `kCells=10` 有跨细胞平滑，不能在本项目“不补 S”的设定下不经说明照搬默认。若未来直接调用，可将平滑作为单独模式审查，或明确采用 `kCells=1`；目前未运行该包或宣称已验证该设置。其 GPL-3 源码仅作为参考，没有移植进当前 Python package。

## scVelo 可借鉴的模块

[官方源码](https://github.com/theislab/scvelo)，本轮审查 commit `f63c0e70596ced2f1bee8cf07e8ab66037cf86b2`；API 检查使用本地 0.3.4。

- `tools/_em_model_core.py` 的 `recover_dynamics` 与 `latent_time`：基因动力学拟合及共享时间。
- `tools/velocity.py`：动力学/稳态响应及基因选择。适配器明确使用 dynamical likelihood≥0.1，并关闭附加稳态 R² 筛选，避免把稳态假设重新加入动力学模式。
- `tools/velocity_graph.py`：速度与候选状态差的方向一致性，用于转移图。
- `preprocessing/moments.py`：邻域一阶/二阶矩。当前适配器不调用此函数，使用 `use_raw=True`，不替换分子观测。

参考官方 [recover_dynamics](https://scvelo.readthedocs.io/en/stable/scvelo.tl.recover_dynamics.html)、[latent_time](https://scvelo.readthedocs.io/en/latest/scvelo.tl.latent_time.html) 接口。实现按实际基因层 precursor/mature counts 检查，不把已剪接异构体的分支比例重新命名为 U/S。

## IsoVelo 的核查边界

[ASHG 2025 会议摘要](https://www.ashg.org/wp-content/uploads/2025/10/2025-Poster-Abstracts.pdf)包含 *Deep Generative Modeling of Isoform Aware RNA Velocity in Long-read Single-cell RNA-seq*，介绍 IsoVelo。检索线索曾指向 `shizhuom/IsoVelo`，本轮 GitHub 返回 404，clone 无法匿名获取。因而尚未审查它的可运行实现，不作为已经可接入的依赖；该访问结果也不能证明工具不存在或永久不可用。

## 实现决策

保留两类独立、可比较的校准：

1. **成熟剪接组成校准**：直接使用事件分支分子计数与有效分母，辅助 G 状态排序；不是 RNA velocity。
2. **动力学定向校准**：只有真正的 precursor/mature 输入及量化来源通过检查时，调用 `splicefate.velocity`，导出所有事件共享的 `latent_time`，再作为可选输入。失败保持空值和 QC 标记。

Region Velocity 提供后续 long-read 区域动力学的实现依据，但需独立 region quantification adapter，不把已有 AS/IR 矩阵冒充区域统计。当前不把当前数据条件强加给通用模型，也不宣称 scVelo 或 Region Velocity 能自动推断每个内含子的真实切除顺序。

## 验证范围

历史版本（当时含基因块排除，当前已移除）的适配器对独立生成的 100 细胞、30 基因 U/S 动力学 Poisson 数据进行了 API 检查，两组基因排除各有 16/14 个动力学支持基因，均完成 velocity graph 和 latent_time。结果用于确认代码路径可执行，不是长读长生物学准确率验证。原始合成输入、日志和输出位于 `/tmp/splicefate_velocity_synthetic`，持久摘要记录在 `runs/implementation_validation/scvelo_adapter_check.json`。
