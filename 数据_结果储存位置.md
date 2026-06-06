# MDD Genetic Risk Analysis Using scDRS and GSMap

## Project Structure

### scDRS 分析

#### 1. 输入数据

Annotated single-cell datasets:

```bash
#male
/biostack/home/suchenyi/scdrs/raw_hg_combined.h5ad
#female
/biostack/home/suchenyi/graduate/scdrs/Annotated_Female1.h5ad
```

GWAS summary statistics:

```bash
/biostack/home/suchenyi/graduate/gwas/scdrs_hg38/
```

#### 2. 分析结果

scDRS step1-compute_scores:

```bash
#male
/biostack/home/suchenyi/graduate/scdrs/score_hg38/
#female
/biostack/home/suchenyi/graduate/scdrs/female_score/
```

scDRS step2-downstream:

```bash
#male
/biostack/home/suchenyi/graduate/scdrs/downstream_hg38/
#female
/biostack/home/suchenyi/graduate/scdrs/female_downstream/
```


#### 3.可视化结果
结果目录：
```text
/biostack/home/suchenyi/graduate/scdrs/results/
```

#### （1）跨性状遗传关联分析
```text
1.dual_panel_trait_correlation.pdf
```
基于 scDRS 与 GSMap 构建的性状相关性热图，用于比较精神疾病、行为性状及代谢性状之间的遗传关联模式。

---
#### （2）脑区特异性遗传富集分析
```text
/biostack/home/suchenyi/graduate/scdrs/results/2.region/
```
包含不同脑区的 MCP 热图、ΔMCP 热图以及 BA46 脑区重点分析结果，用于评估 MDD 相关遗传风险在脑区间的分布差异。
主要文件：
```text
2.1scDRS_deltaMCP_conSection_lambdaGC_heatmap
2.2assoc_mcz_region_alltraits_BA46_highlighted
scDRS_triangle_MCP_heatmap_lambdaGC_section
```
---

#### （3）MDD性别差异分析
```text
/biostack/home/suchenyi/graduate/scdrs/results/3.male vs female/
```
比较男性与女性样本中的遗传风险富集模式，包含 MCP 热图、ΔMCP 热图及性别差异热图。
主要文件：
```text
scDRS_deltaMCP_female_male_heatmap
scDRS_deltaMCP_sex_diff_heatmap
scDRS_MCP_male_female_sidebyside
```
---

#### （4）女性特异性分析
```text
/biostack/home/suchenyi/graduate/scdrs/results/4.1female/
```
针对女性数据开展细胞类型、细胞亚型及脑区水平的遗传风险富集分析。
包括：
* Cell Type（CT）
* Cell Subtype（CTS）
* Brain Region（Section）
三个层次的 MCP 及 ΔMCP 结果。

---

#### （5）男性特异性分析
```text
/biostack/home/suchenyi/graduate/scdrs/results/4.2male/
```
针对男性数据开展细胞类型、细胞亚型及脑区水平的遗传风险富集分析。

---

#### （6）跨性状网络分析

```text
MDD_chord_v3.pdf
Psych_Metabolic_Sankey.html
Sankey_C.html
Sankey_S.html
```
用于展示 MDD 与其他精神疾病、行为性状及代谢性状之间的共享遗传基础。
包括：
* Chord Diagram（弦图）
* Sankey Diagram（桑基图）

---

#### （7）单细胞遗传风险空间分布
```text
/biostack/home/suchenyi/graduate/scdrs/umap_hg38/
```
包含全部 GWAS 性状的 UMAP 可视化结果。
每个 PDF 对应一个性状在单细胞图谱中的遗传风险富集分布，例如：
```text
64human_MDD_deepblue_bar.pdf
64human_SCZ_deepblue_bar.pdf
64human_Neuroticism_deepblue_bar.pdf
64human_E-Smoking_deepblue_bar.pdf
64human_AD_deepblue_bar.pdf
```
用于观察疾病相关细胞群体及其空间分布特征。

---


### GSMap分析

#### 1.输入数据

原始空转数据:
```bash
#human
/biostack/home/suchenyi/graduate/gsmap/raw_data/
#mouse
/biostack/home/suchenyi/gsmap/st/
```

GWAS summary statistics:
```bash
/biostack/home/suchenyi/gwas/gsmap_hg38/
```
gsMap工作及结果输出目录：
```bash
#human
/biostack/home/suchenyi/graduate/gsmap/human_h5ad
#mouse
/biostack/home/suchenyi/gsmap/
```


#### 2. GSMap可视化结果

主要结果目录：

```text
/biostack/home/suchenyi/graduate/gsmap/human_results/
/biostack/home/suchenyi/graduate/gsmap/mouse_results/
```
#### （1）人脑空间遗传风险分布
基于 BA12 和 BA46 脑区空间转录组数据，展示 MDD 遗传风险在不同脑区及细胞类型中的空间分布。
主要结果包括：

```text
gsmap_spatial_BA46_rep1
gsmap_spatial_by_rep
gsmap_celltype_spatial
gsmap_2x2_spatial
```
以及 BA12、BA46 不同条件下的风险评分分布小提琴图。

#### （2）小鼠慢性应激模型空间分析
```text
/biostack/home/suchenyi/graduate/gsmap/mouse_results/gsmap_spatial_plots/
```
分析 mPFC、HIP、VTA、DRN 和 Stri 等脑区在 Control、3D、7D、Sus 和 Res 状态下的遗传风险空间分布变化。

#### （3）细胞类型富集分析
```text
/biostack/home/suchenyi/graduate/gsmap/mouse_results/celltype_violin/
/biostack/home/suchenyi/graduate/gsmap/mouse_results/OR/
```
包括：
* 细胞类型风险评分小提琴图
* OR（Odds Ratio）富集分析
* BA12 与 BA46 细胞类型比较
用于识别 MDD 相关风险富集的关键细胞群体。

#### （4）跨性状遗传结构分析
```text
/biostack/home/suchenyi/graduate/gsmap/mouse_results/traits_corr/
```
包括：
* Trait Correlation Heatmap
* Hierarchical Clustering
* PCA
用于比较不同疾病和行为性状之间的共享遗传结构。

#### （5）应激轨迹分析

```text
/biostack/home/suchenyi/graduate/gsmap/mouse_results/trajectory/
```
基于 Control → 3D → 7D → Sus/Res 的时间轨迹，分析 MDD、SCZ、Neuroticism、E-Smoking、IQ 等性状遗传风险在不同细胞类型中的动态变化。

---
