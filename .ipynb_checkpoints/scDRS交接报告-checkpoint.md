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

*   `1.submit_scdrs1.sh`: Shell脚本，概述了主要的运行流程，包含了GWAS数据处理的关键步骤和简单的scDRS运行示例。可作为理解整体流程的参考。
*   `2.submit_plot.sh`
*   `2.1.submit_scplot.sh`
*   ``
*   ``
*   `scdrs_run1.ipynb`: **核心文件**，Jupyter Notebook格式，包含了完整的scDRS运行代码、参数设置、以及结果的可视化绘图代码。**主要的分析和绘图都在此文件中完成。**
*   `umap_sc.py`: Python补充脚本，用于对特定的单细胞数据（如Zeisel 2015）进行预处理，可作为准备scRNA-seq输入数据 `.h5ad` 格式的参考。
*   `umap_plot.py`：
*   `4.1.lambdaGC.py`：
*   `lambdaGC_barplot.py`：
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

完整的分析代码、参数调整和可视化均在 `scdrs_run1.ipynb` Jupyter Notebook中。

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

    * 热图展示了 细胞类型 [dominant_labels,dominant_labels_short] 与 不同性状 的平均 scDRS MCP 关联值。

    * 每个格子表示某细胞类型与某性状的关联程度，颜色越深表示关联越强（0–1 归一化）。

    * 聚类分析：

        * 行聚类 (row cluster)：将细胞类型按照它们在不同性状上的关联模式相似性进行分组。

        * 列聚类 (col cluster)：将性状按照它们在不同细胞类型上的关联模式相似性进行分组。
     
5. **scDRS 群体水平热图（带 λGC）**

    * λGC 柱状图

        * 顶部的灰色/红色柱状图显示每个性状的 λGC 值（基因组学中用于衡量系统性偏倚或多重检验膨胀的指标）。

        * λGC 越高，说明对应的 GWAS 性状统计量可能存在轻微膨胀，需要注意。

        * 颜色深浅对应 λGC 大小，方便直观比较不同性状的偏倚程度。



