# scDRS交接报告

scDRS是一款用于联合分析单细胞RNA测序（scRNA-seq）数据或空间转录组数据与GWAS（全基因组关联研究）疾病相关性状数据的软件。

主要的分析步骤包括：
1.  手动收集GWAS summary statistics数据，并进行**多步格式转换**以满足scDRS和中间工具（MAGMA）的要求。
2.  执行scDRS分析流程，主要包括两个核心部分：
    *   计算每个细胞相对于GWAS性状的疾病相关性得分（disease relevance score）。
    *   评估在特定细胞类型内部，这种疾病相关性得分的异质性（heterogeneity）。

关键点：前期准备GWAS数据阶段需要手动收集且格式转换步骤繁琐，比较耗时。后续的scDRS计算流程相对较快（使用Jupyter Notebook执行），对示例中的13个性状，预计半天内可完成计算和初步可视化。

参考链接：
*   官网：[https://github.com/martinjzhang/scDRS](https://github.com/martinjzhang/scDRS)
*   文档：[Quick start — scDRS 1.0.3 documentation](https://martinjzhang.github.io/scDRS/notebooks/quickstart.html)
*   中文流程1（官网中文版）：[Post-GWAS： single-cell disease relevance score (scDRS) 分析 - 知乎](https://zhuanlan.zhihu.com/p/592128325)
*   中文流程1-2（预处理 - MAGMA）：[基于 MAGMA 的 gene-based 关联分析研究 - 橙子牛奶糖 - 博客园](https://www.cnblogs.com/chenwenyan/p/14628970.html)
*   中文流程2（完整代码）：[scDRS 单细胞RNA与GWAS联合分析（全网最详细！！！）](https://mp.weixin.qq.com/s/iZsbRgyjI1T5q3moEbTvQA)
*   中文流程3（完整代码）：[进行单细胞和GWAS的scDRS分析（三）：进行scDRS主体分析](https://mp.weixin.qq.com/s/I52r0BHkuDt1F9h_IVnurw)
*   中文流程3-2（预处理 - MAGMA）：[使用python进行单细胞和GWAS的scDRS分析（二）——MAGMA工具使用](https://mp.weixin.qq.com/s/mWcucUSt6-PxcrAURUd8-A)

补充材料：scRNA-GWAS-software-Collect.xlsx，summaryStats.pdf（网上GWAS数据库pdf）

---

## 详细分析流程

### 整体流程代码文件概览

以下是与scDRS分析相关的主要脚本和文件列表：

*   `1.submit_scdrs1.sh`: 主要的运行流程一scdrs compute-score。
*   `2.submit_plot.sh`：整体结果画图批量运行代码
*   `2.1.submit_scplot.sh`：分疾病组与控制组结果画图批量运行代码
*   `3.scdrs_downstream.sh`：主要的运行流程二整体下游分析
*   `scdrs_run1-mdd3.ipynb`: Jupyter Notebook格式，包含了（单性状）完整的scDRS运行代码、参数设置、以及部分结果的可视化绘图代码。
*   `01_scDRS_male.ipynb`: **核心文件**，Jupyter Notebook格式，包含了（多性状）男性结果的可视化绘图代码。**主要的scDRS结果绘图都在此文件中完成。**
*   `01_scDRS_female.ipynb`: **核心文件**，Jupyter Notebook格式，包含了（多性状）女性结果的可视化绘图代码。**主要的男女性比较绘图都在此文件中完成。**
*   `umap_sc.py`: 分疾病组与控制组结果画图代码
*   `umap_plot.py`：整体结果画图代码
*   `4.1.lambdaGC.py`：多性状的λGC计算代码
*   `lambdaGC_barplot.py`：对多性状的λGC结果可视化脚本
*   `tsv2snpp.R`: R脚本，用于将初步整理好的GWAS summary文件（`.tsv`格式）转换为MAGMA软件所需的SNP P值文件（`.snpp`格式），并提取样本量信息。
*   `generate_magma1.sh`: Shell脚本，用于批量生成运行MAGMA所需的命令，将`.snpp`文件和参考基因组信息结合，计算基因水平的关联P值（生成`.genes.out`文件）。
*   `genout2gs.R`: (在此报告中提及，用于后续步骤) R脚本，用于将多个MAGMA输出的`.genes.out`文件合并，并整理成scDRS格式化工具所需的输入格式 (`geneset.pval`)。

### 需手动准备和格式化的GWAS数据

**重要提示：** 此步骤与gsMap的GWAS数据准备类似，但后续格式转换更为复杂，需要经过 **GWAS -> .tsv -> .snpp -> .genes.out -> .gs** 的多步流程。

1.  **数据来源与选择**：
    *   同gsMap报告所述，主要从 [GWAS Catalog](https://www.ebi.ac.uk/gwas/) 或原始文献补充材料获取Full Summary Statistics。
    *   选择感兴趣的性状进行下载。

2.  **格式化转换流程【核心难点】**：

    *   **步骤 2.1: 下载GWAS -> 标准化 `.tsv` 文件**
        *   与gsMap类似，使用`awk`等工具处理下载的原始文件，生成一个包含关键列的 `.tsv` 文件。
        *   **关键列：** 对于后续的MAGMA和scDRS流程，最重要的列是 `SNPID` (或`RSID`)、`P` (P值) 和 `N` (样本量)。其他列如 `CHR`, `POS`, `A1`, `A2` 等在初步整理时也建议保留，但 `BETA`, `SE` 对此流程不是必需的。
        *   示例 `.tsv` 文件 (`trait_Depression.tsv`)：
    
    ```
    SNPID   CHR     POS     RSID    A1      A2      OR      BETA    SE      P       N
    rs12049979      12      110352987       rs12049979      A       C       0.9998  -0.00020002     0.0061  0.977   294323
    rs1807335       12      113545165       rs1807335       T       C       1.0026  0.00259663      0.0044  0.5467  294322
    rs921063        12      118274349       rs921063        T       C       1.00642 0.00639948      0.0036  0.07328 294323
    ```

*   **步骤 2.2: `.tsv` -> `.snpp` (使用 `tsv2snpp.R`)**
    
    *   目的：转换为MAGMA软件可以识别的SNP P值输入格式。
    *   运行 `tsv2snpp.R` 脚本。该脚本会：
        *   读取目录中所有的 `.tsv` 文件。
        *   提取 `RSID` 和 `P` 列。
        *   为每个 `.tsv` 文件输出一个对应的 `.snpp` 文件。
        *   同时，计算每个 `.tsv` 文件中 `N` 列的众数（因为N可能略有不同），并输出到一个汇总文件 `summary_Nsample.tsv` 中，供下一步使用。
    *   示例 `summary_Nsample.tsv` 文件：
        ```
        file             N_mode
        trait_A-ADHD.tsv  58286
        trait_BMI.tsv     795640
        trait_BP.tsv      413466
    ```
    
*   **步骤 2.3: `.snpp` -> `.genes.out` (使用 MAGMA，通过 `generate_magma1.sh` 批量执行)**
    *   目的：使用MAGMA软件，基于SNP P值 (`.snpp`) 和参考基因组注释 (`g1000_eur_qc.genes.annot`)，计算基因水平的关联P值。
    *   运行 `bash generate_magma1.sh`。该脚本会：
        *   读取 `summary_Nsample.tsv` 获取每个性状对应的样本量 `N`。
        *   为每个 `.snpp` 文件生成一条MAGMA命令，指定输入 `.snpp` 文件、样本量 `N`、参考基因组文件 (`--bfile g1000_eur_qc`)、基因注释文件 (`--gene-annot`) 和输出文件前缀。
        *   将所有生成的MAGMA命令写入 `magma1.sh` 文件。
    *   `generate_magma1.sh` 脚本内容：
        ```bash
        #!/bin/bash
        # generate_magma1.sh: Creates a script (magma1.sh) containing MAGMA commands for gene-level association analysis.
        #################################
        output_file="magma1.sh"
        > "$output_file" # Clear or create the output script file

        # Read summary_Nsample.tsv to create a mapping from trait prefix to sample size (N)
        declare -A prefix_to_N
        while IFS=$'\t' read -r file N; do
            if [[ "$file" == *.tsv && "$file" != "file" ]]; then # Skip header
                prefix=${file%.tsv}  # Remove .tsv suffix
                prefix_to_N[$prefix]=$N
            fi
        done < summary_Nsample.tsv

        # Iterate through all .snpp files and generate MAGMA commands
        for snpp_file in *.snpp; do
            if [[ -f "$snpp_file" ]]; then
                prefix=${snpp_file%.snpp}  # Remove .snpp suffix
                if [[ -n "${prefix_to_N[$prefix]}" ]]; then
                    N=${prefix_to_N[$prefix]}
                    # Construct the MAGMA command
                    cmd="magma --bfile g1000_eur_qc --pval \"$snpp_file\" N=$N --gene-annot g1000_eur_qc.genes.annot --out \"$prefix\""
                    echo "$cmd" >> "$output_file"
                else
                    echo "Warning: No N value found for prefix $prefix in summary_Nsample.tsv" >&2
                fi
            fi
        done
        echo "MAGMA commands written to $output_file. Run 'bash $output_file' or 'parallel -j <num_jobs> < $output_file' to execute."
        ```
    *   **执行MAGMA命令**: 由于MAGMA计算相对独立且占用资源不大，可以使用 `parallel` 工具并行执行 `magma1.sh` 中的命令：
        ```bash
        parallel -j 13 < magma1.sh # 示例：使用13个核心并行运行
        ```
        这将为每个性状生成一个 `.genes.out` 文件。

*   **步骤 2.4: `.genes.out` -> `.gs` (使用 `genout2gs.R` 和 `scdrs munge-gs`)**
    *   目的：将MAGMA生成的基因水平P值 (`.genes.out`) 转换为scDRS可以直接使用的基因集文件 (`.gs`)。
    *   **第一步：合并与整理 (使用 `genout2gs.R`)**:
        *   运行 `Rscript genout2gs.R`。这个脚本的作用是：
            *   读取所有生成的 `.genes.out` 文件。
            *   提取每个基因在每个性状下的P值。
            *   将结果合并成一个单一的、基因×性状的P值矩阵文件，命名为 `geneset.pval`。
        *   示例 `geneset.pval` 文件 (`head` 输出)：
            ```
            GENE    A-ADHD  BMI     BP      CoPe    CRP     Depression      Depression2     E-Smoking       Insomnia        IQ      MDD     Neuroticism     SCZ
            A1BG    0.97899 0.00022439      0.91978 0.22365 0.43348 0.34535 0.58203 0.29819 0.15962 0.09873 0.061758        0.073231        0.60464
            A1CF    0.59911 0.71002 0.0056546       0.041815        0.030248        0.28988 0.95617 0.42198 0.18979 0.050854        0.67789 0.35455 0.55488
            ```
    *   **第二步：最终格式化 (使用 `scdrs munge-gs`)**:
        *   使用scDRS自带的 `munge-gs` 命令完成最后转换。
        ```bash
        # Convert the combined p-value file to the final .gs format required by scDRS
        scdrs munge-gs --out-file "geneset1.gs" --pval_file "geneset1.pval" --weight zscore --n-max 1000
        ```
        *   `--weight zscore`: 指定使用Z-score权重。
        *   `--n-max 1000`: 可能与基因集大小或选择有关（具体参考scDRS文档）。
        *   最终得到 `geneset1.gs` 文件，作为scDRS分析的输入。

### 主要的scDRS分析流程 (在 `scdrs_run1.ipynb` 中执行)

完整的分析代码、参数调整和部分可视化均在 `scdrs_run1-mdd3.ipynb` Jupyter Notebook中。

**重要前提：** 输入的单细胞/空间转录组表达矩阵 (`.h5ad` 文件中的表达层) **不能包含负数值**。如果数据经过了可能产生负数的标准化（如某些 `scanpy.pp.scale` 操作），需要使用原始计数或经过对数转换（如 `log1p`）且非负的数据层。

**分析结果可视化示例：**

1.  **scDRS 核心结果：细胞水平疾病相关性得分**

    * 展示了每个细胞（点）根据其基因表达谱计算得到的与特定GWAS性状（如BP - 血压）的关联得分（scDRS score），并按细胞类型（annotation）着色。得分越高表示该细胞的表达模式与该性状的遗传风险关联越强。

2.  **群体水平关联与异质性总结**
    *   **图例解释：**
        *   **热图颜色 (Heatmap color)**：表示在某个细胞类型-性状对中，表现出显著关联（即scDRS得分显著非零）的细胞所占的比例。颜色越深，比例越高。
        *   **方块 (Squares)**：标记显著的细胞类型-性状关联（在所有细胞类型和性状中进行多重检验校正后，FDR < 0.05）。表示该细胞类型整体上与该性状显著相关。
        *   **叉号 (Cross symbols)**：标记在某个特定细胞类型内部，不同细胞间的scDRS得分存在显著的异质性。表示即使该细胞类型整体与性状相关，其内部细胞的关联程度也有显著差异。
        *   注：后面由于异质性结果不明显，故下一步画图把异质性去掉了，只关注关联性
3. **热图展示scDRS MCP 关联值**

    * trait × condition
    * trait × dominant_labels
    * trait × dominant_labels_short
    * section × dominant_labels
    * section × dominant_labels_short
    * condition × dominant_labels_short
    * condition × dominant_labels
为了后面画BA×celltype的热图，需要把BA×celltype的每个组合作为下游分析的分类依据进行下游分析得到每个组合的MCP值

```bash
scdrs perform-downstream \
    --h5ad-file adata_BA_celltype.h5ad \
    --score-file mdd3_analysis/mdd.full_score.gz \
    --group-analysis "BA_celltype" \
    --flag-filter-data False \
    --flag-raw-count True \
    --out-folder mdd3_analysis/
```
      
输入的h5ad文件只是在原来的raw_hg.h5ad上新建了一列BA_celltype，其实就是每行的section和dominant_labels(或dominant_labels_short)组合
还有condition和dominant_labels(或dominant_labels_short)组合，一共四组

4. **scDRS 群体水平关联热图**
* **`01_scDRS_male.ipynb`cell-6+7+8** *
* **`01_scDRS_female.ipynb`cell-2+3** *
    * 热图展示了 细胞类型 [dominant_labels,dominant_labels_short] 与 不同性状 的平均 scDRS MCP 关联值。

    * 每个格子表示某细胞类型与某性状的关联程度，颜色越深表示关联越强（0–1 归一化）。

    * 聚类分析：

        * 行聚类 (row cluster)：将细胞类型按照它们在不同性状上的关联模式相似性进行分组。

        * 列聚类 (col cluster)：将性状按照它们在不同细胞类型上的关联模式相似性进行分组。

    * λGC 柱状图

        * 顶部的灰色/红色柱状图显示每个性状的 λGC 值（基因组学中用于衡量系统性偏倚或多重检验膨胀的指标）。

        * λGC 越高，说明对应的 GWAS 性状统计量可能存在轻微膨胀，需要注意。

        * 颜色深浅对应 λGC 大小，方便直观比较不同性状的偏倚程度。
     
5. **scDRS MCP 上下三角差异关联热图**
* **`01_scDRS_male.ipynb`cell-13+14+15** *
* **`01_scDRS_female.ipynb`cell-6+7+8+9** *
   * 热图展示了不同细胞类型/脑区（dominant_labels 或 dominant_labels_short 或 region）在 MDD（S）与 Control（C）条件下的 scDRS 群体水平关联差异。
   * 每个格子由两个三角形组成：
      * 上三角表示 MDD 组（S）的 assoc.MCP；
      * 下三角表示 Control 组（C）的 assoc.MCP。

    * 颜色表示关联显著性水平（Assoc. MCP）：
        * 红色表示 MCP 值较低，即遗传风险与该细胞类型的关联更显著；
        * 蓝色表示 MCP 值较高，即关联较弱；
        * 颜色梯度反映不同条件下细胞类型与性状之间的关联强度差异。
6.  **scDRS ΔMCP 状态差异关联热图**
* **`01_scDRS_male.ipynb`cell-11+10+9** *
* **`01_scDRS_female.ipynb`cell-4+5+11** *
    * 热图展示了不同细胞类型/脑区（dominant_labels 或 dominant_labels_short 或 region）在 MDD（S）与 Control（C）条件下的 scDRS 群体水平关联差异。
    * 每个格子表示某一细胞类型与某一 GWAS 性状的 ΔMCP 值。

    * 颜色表示疾病状态相对于对照状态的关联变化方向和幅度：
        * 红色表示 ΔMCP > 0，即 MDD 组 MCP 高于 Control 组；
        * 蓝色表示 ΔMCP < 0，即 MDD 组 MCP 低于 Control 组；
        * 颜色越深表示两组之间的差异越大。
    * 由于 MCP 值越小代表遗传风险富集越显著，因此：
        * ΔMCP < 0 表示 MDD 组的遗传风险关联增强；
        * ΔMCP > 0 表示 MDD 组的遗传风险关联减弱。
7. **scDRS 与 GSMap 性状相关性双面板热图**
* **`01_scDRS_male.ipynb`cell-16** *
    * 热图展示了不同 GWAS 性状之间遗传风险富集模式的相似性，并比较 scDRS 与 GSMap 两种方法所得结果的一致性。
    * 图中包含两个面板：
        * **Panel A（scDRS）**：基于各性状在不同细胞类型中的 assoc.MCP 值计算性状间相关性；
        * **Panel B（GSMap）**：基于各性状在空间转录组位点中的 z_mean 值计算性状间相关性。

    * 每个格子表示两种性状之间的 Pearson 相关系数（Pearson's r）。
    * 颜色表示性状间关联模式的相似程度：
        * 红色表示正相关，即两种性状在细胞类型或空间位点中的遗传风险富集模式相似；
        * 蓝色表示负相关，即两种性状呈现相反的富集模式；
        * 颜色越深表示相关性绝对值越高。
    * 聚类分析：
        * 性状排序基于 Panel A（scDRS）的层次聚类结果。
        * 聚类依据为各性状在所有细胞类型中的关联谱（assoc.MCP），将具有相似细胞类型富集模式的性状聚集在一起。
        * 左侧树状图（dendrogram）展示了性状间的层次聚类关系，不同性状簇反映其潜在共享的遗传风险特征。
    * 双面板比较：
        * 两个面板采用完全相同的性状排序，便于直接比较 scDRS 与 GSMap 所识别的性状相关结构。
        * 若两个面板呈现相似的相关性模块，说明细胞类型层面与空间层面的遗传风险富集模式具有较好一致性。
        * 若某些性状簇仅在其中一个面板中表现出较强相关性，则提示其细胞类型特异性与空间分布特征可能存在差异。
    * 颜色条（Color Bar）：
        * 右侧颜色条表示 Pearson 相关系数范围。
        * Pearson's r 取值范围为 −1 至 1：
            * r > 0 表示正相关；
            * r < 0 表示负相关；
            * |r| 越接近 1 表示相关性越强；
            * r ≈ 0 表示两种性状间缺乏明显相关关系。
    * 为提高结果可靠性，仅纳入 λGC ≤ 2 的 GWAS 性状进行分析，以减少统计量膨胀对相关性结构的潜在影响。
8. **scDRS ΔMCP 性别差异关联热图**
* **`01_scDRS_female.ipynb`cell-14** *
    * 热图展示了不同细胞类型（dominant_labels_short）在 female与male条件下的 scDRS 群体水平关联差异。
    * 每个格子表示某一细胞类型与某一 GWAS 性状的 ΔMCP (female-male)值。

    * 颜色表示疾病状态相对于对照状态的关联变化方向和幅度：
        * 红色表示 ΔMCP > 0，即女性组 MCP 高于男性组；
        * 蓝色表示 ΔMCP < 0，即女性组 MCP 低于男性组；
        * 颜色越深表示两组之间的差异越大。
    * 由于 MCP 值越小代表遗传风险富集越显著，因此：
        * ΔMCP < 0 表示 女性组的遗传风险关联更强；
        * ΔMCP > 0 表示 男性组的遗传风险关联更强。

9. **scDRS MCP 上下三角性别差异关联热图**
* **`01_scDRS_female.ipynb`cell-12** *
   * 热图展示了不同细胞类型（dominant_labels）在 female与male条件下的 scDRS 群体水平关联差异。
   * 每个格子由两个三角形组成：
      * 上三角表示男性的 assoc.MCP；
      * 下三角表示女性的 assoc.MCP。

10. **脑区-性状关联 z-score 分布图**
* **`01_scDRS_male.ipynb`cell-12** *
    * 图中展示了不同脑区在 Control（C）与 MDD（S）条件下的遗传风险关联强度分布。
    * 关联强度采用 association z-score（assoc_mcz）表示，用于衡量特定脑区与 GWAS 性状之间的遗传风险富集程度。

    * 箱线图（Boxplot）：
        * 每个脑区包含两组箱线图，分别对应 Control（蓝色）和 MDD（红色）条件。
        * 箱体表示四分位距（IQR），中线表示中位数，反映该脑区所有性状关联强度的总体分布特征。
        * 箱线图用于比较不同条件下脑区整体遗传风险富集水平的变化趋势。
    * 散点图（Trait-level points）：
        * 每个散点代表一个 GWAS 性状在对应脑区中的 association z-score。
        * 彩色散点表示重点关注的神经精神相关性状，包括 SCZ、MDD、Neuroticism、Insomnia、IQ 等；
        * 灰色散点表示其余 GWAS 性状。
        * 散点分布能够反映不同性状之间的异质性及其在各脑区中的遗传风险富集差异。
    * 中位数变化连线：
        * 每个脑区均绘制从 Control 到 MDD 的中位数连线，用于展示整体关联强度的变化方向。
        * 连线起点表示 Control 条件下所有性状 assoc_mcz 的中位数；
        * 连线终点表示 MDD 条件下所有性状 assoc_mcz 的中位数；
        * 连线向上表示疾病状态下整体关联增强，向下表示整体关联减弱。
    * BA46 脑区高亮显示：
        * BA46 脑区采用红色连线和标记进行突出显示。
        * 该脑区在疾病与对照条件之间表现出较明显的关联强度变化，因此作为重点研究区域进行展示。
        * 连线旁标注的数值表示对应条件下 assoc_mcz 的中位数。
    * 参考线：
        * 水平虚线表示 z-score = 0。
        * z-score > 0 表示正向富集；
        * z-score < 0 表示负向富集；
        * 绝对值越大表示关联强度越高。






