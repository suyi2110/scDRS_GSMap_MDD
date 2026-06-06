# gsMap交接报告

gsMap是用于联合分析空间转录组（Spatial Transcriptomics）数据和GWAS（全基因组关联研究）疾病相关性状数据的软件。

主要的分析步骤包括：
1.  手动收集和格式化GWAS summary statistics数据（根据目标性状选择）。
2.  使用gsMap对空间转录组数据进行预处理（例如切片取平均）。
3.  运行gsMap核心分析流程（分步执行：寻找潜在表征 -> 计算GSS -> 生成LD score -> 空间LDSC分析）。

关键点：前期准备GWAS数据阶段需要手动收集，比较耗时。同时gsMap运行时需要的计算资源较多（可能会占用大量CPU），建议采用分步流程（step-by-step），便于中断后从失败步骤继续运行。

参考链接：
*   官网：[github: gsMap](https://github.com/JianYang-Lab/gsMap)
*   文档：[gsMap documentation](https://yanglab.westlake.edu.cn/gsmap/document/software)

补充材料：gsMap结果报告解读v2.docx，GWAS-Sumstats-collect_20250409.xlsx，

---

## 详细分析流程

### 整体流程代码文件概览

以下是与gsMap分析相关的主要脚本文件列表：
*   `submit_gsmap_tasks_cpu2.sh`: 用于在后台提交批量任务的辅助脚本。
*   `submit_spatial_ldsc_array.sh`: 用于通过`sbatch`（SLURM调度系统）提交批量任务的辅助脚本。
*   `01-slice-mean-run-hg.sh`: **步骤1** - 对人类空间转录组切片进行重复样本取平均。
*   `02-find-latent-run.sh`: **步骤2** - 生成寻找潜在表征（latent representation）的任务列表。
*   `03-latent2gene-run-hg.sh`: **步骤3** - 生成计算GSS（Gene Specific Score）的任务列表（人类）。
*   `04-generate-ldscore.sh`: **步骤4** - 生成LD score的任务列表。
*   `05-spatial-ld&domain.sh`: **步骤5 & 6** - 生成空间LDSC（Spatial LD Score Regression）分析和报告生成的任务列表。

### 需手动准备的GWAS数据

gsMap分析需要格式正确的GWAS summary statistics数据。以下是准备步骤：

1.  **数据来源**：
    
    *   主要来源是 [GWAS Catalog](https://www.ebi.ac.uk/gwas/)。
    *   如果Catalog不提供Full Summary Statistics，需要到原始研究论文中查找补充材料，数据可能托管在其他数据库或个人网站（如FigShare）。
    
2.  **数据选择**：
    *   根据研究目标，在GWAS Catalog或其他来源搜索感兴趣的性状（例如"Depression"）。
    *   下载**Full Summary Statistics**文件。这些文件通常较大（几百MB到1GB），包含该研究中所有SNP（包括显著和不显著）的统计信息。每个文件代表一个独立的GWAS分析。
    
    ![](E:\00-SMU-PC\PhD-work\00-report-email\md-pictures\image-20250428094413289.png)
    
    ![image-20250428094459510](E:\00-SMU-PC\PhD-work\00-report-email\md-pictures\image-20250428094459510.png)
    
    ![image-20250428094704885](E:\00-SMU-PC\PhD-work\00-report-email\md-pictures\image-20250428094704885.png)
    
    注意下载的时候要选择较大的文件，有些并不包括所有的SNP。
    
    ![image-20250428094806965](E:\00-SMU-PC\PhD-work\00-report-email\md-pictures\image-20250428094806965.png)
    
    ![image-20250428094855306](E:\00-SMU-PC\PhD-work\00-report-email\md-pictures\image-20250428094855306.png)
    
3. **关键列确认**：
   *   确保下载的数据包含以下关键列（或可以通过已有列转换得到）。列名可能不同，需要手动确认对应关系：
       *   `SNPID` / `RSID`: SNP标识符
       *   `CHR`: 染色体
       *   `POS` / `BP`: 物理位置
       *   `A1`:效应等位基因 (Effect Allele)
       *   `A2`:非效应等位基因 (Other Allele)
       *   `PVAL` / `P`: P值
       *   `BETA` / `OR`: 效应值 (通常用BETA，没有BETA时有时可用log(OR)代替)
       *   `SE`: 标准误
       *   `N`: 样本量 (非常重要，每个GWAS数据不同)
   *   **格式转换参考**：[GWAS summary statistics 文件格式处理 - AllenW - 博客园](https://www.cnblogs.com/chenwenyan/p/14481976.html)

4. **格式化 (两步)**：
   *   **第一步：初步整理**：使用`awk`或其他工具，将下载的原始文件整理成包含上述关键列、以制表符分隔的`.tsv`文件，并统一列名。
       *   **示例** (`GCST90019499_buildGRCh37.tsv` 处理为 `CRP.tsv`)：
         ```bash
         # 注意：这里的列索引($1, $7等)和N值需要根据你下载的具体文件进行调整
         awk 'BEGIN {OFS="\t"}
              NR == 1 {
                  # 输出gsMap期望的标准表头
                  print "SNPID", "CHR", "POS", "RSID", "A1", "A2", "P", "BETA", "SE", "N";
                  next
              }
              {
                  # 根据原始文件的列定义，提取对应信息
                  A1 = $1      # effect_allele
                  A2 = $7      # other_allele
                  POS = $8     # base_pair_location
                  BETA = $9    # beta
                  CHR = $10    # chromosome
                  P = $11      # p_value
                  SE = $12     # standard_error
                  RSID = $13   # variant_id（通常等同于 SNPID）

                  # 输出整理后的行，N值根据文献或数据描述填写（这里示例为固定值）
                  print RSID, CHR, POS, RSID, A1, A2, P, BETA, SE, 355127
              }' GCST90019499_buildGRCh37.tsv > CRP.tsv
         ```
   *   **第二步：gsMap内部格式化**：使用gsMap自带的工具将`.tsv`文件转换为其内部使用的压缩格式 (`.sumstats.gz`)。
       ```bash
       gsmap format_sumstats --sumstats 'CRP.tsv' --out 'CRP'
       ```
       这将生成最终可用于gsMap分析的 `CRP.sumstats.gz` 文件。

### 主要的gsMap分析流程 

**总体策略**：
*   Step 1（切片平均）单独运行。
*   Step 2-4（寻找latent -> 计算GSS -> 生成LD score）合并在一个流程中顺序执行，但每个子步骤内部的任务是并行提交的。
*   Step 5-6（空间LDSC分析）合并运行。

**执行细节**：

1.  **步骤 1: 切片平均 (Slice Mean)**
    *   目的：对来自同一条件、同一区域的重复空间转录组切片数据进行合并取平均，减少数据冗余，得到代表性的表达谱。
    *   脚本：`01-slice-mean-run-hg.sh`
    *   运行命令（后台执行）：
        ```bash
        nohup bash 01-slice-mean-run-hg.sh "/biostack/home/suchenyi/graduate/gsmap/data" "/biostack/home/suchenyi/graduate/gsmap/human_h5ad" > 01-slice-mean-run-hg.log 2>&1 &
        ```
    *   监控：观察 `01-slice-mean-run-hg.log` 文件或使用 `top`/`htop` 检查 `gsmap slice_mean` 进程。

2.  **步骤 2-4: 核心分析流程 (Find Latent -> GSS -> LD Score)**
    *   目的：依次执行gsMap的核心计算步骤，提取空间表达模式（latent representation），计算基因特异性得分（GSS），并为后续与GWAS关联准备LD score。
    *   脚本：依次调用 `02-find-latent-run.sh`, `03-latent2gene-run-hg.sh`, `04-generate-ldscore.sh` 来生成任务列表，并使用 `submit_spatial_ldsc_array.sh` 提交任务
    *   运行命令：
        ```bash
        ################ 步骤 02: 寻找 Latent Representation ################
        bash 02-find-latent-run.sh "/biostack/home/suchenyi/graduate/gsmap/human_h5ad" "/biostack/Data/20250319_CSDS_MDD_annotation/human_stereoseq" "manual_cell_type_reannotated" "run_find_latent_jobs-hg.sh"
        #第 02 步任务数量：
        wc -l < run_find_latent_jobs-hg.sh
        # 使用辅助脚本提交任务，并行度为 10，日志输出到 nohup_logs/02-*.log
        nohup bash submit_gsmap_tasks_cpu2.sh run_find_latent_jobs-hg.sh 10 > nohup.out 2>&1 &

        ################ 步骤 03: 计算 GSS (Gene Specific Score) ################
        bash 03-latent2gene-run-hg.sh "/biostack/home/suchenyi/graduate/gsmap/human_h5ad" "/biostack/Data/20250319_CSDS_MDD_annotation/human_stereoseq" "manual_cell_type_reannotated" "run_latent2gene_jobs-hg.sh"
        nohup bash submit_gsmap_tasks_cpu2.sh run_latent2gene_jobs-hg.sh 10 > nohup.out 2>&1 &
        
        ################ 步骤 04: 生成 LD Score ################
        # 注意：这里需要提供物种对应的 GTF 文件路径
        bash 04-generate-ldscore.sh "/biostack/home/suchenyi/graduate/gsmap/human_h5ad" "/biostack/home/suchenyi/graduate/gsmap/data" "/biostack/home/suchenyi/graduate/gsmap/reference/GRCh38-3.0.0-genes.gtf" "run_generate_ldscore_jobs-hg.sh"
        nohup bash submit_gsmap_tasks_cpu2.sh run_generate_ldscore_jobs-hg.sh 10 > nohup.out 2>&1 &
        ```
    *   **监控与错误处理**：
        *   通过 `nohup_logs/` 目录下的具体日志文件 (`02-*.log`, `03-*.log`, `04-*.log`) 查看每个任务的运行情况。
        *   使用 `pgrep -f gsmap` 或 `htop` 监控 `gsmap` 进程。
        *   如果某个子步骤的任务失败，可以根据 `success_${SPECIES}-*.txt` 文件找到未完成的任务，生成重跑脚本。

3.  **步骤 5-6: 空间 LDSC 分析 (Spatial LDSC & Reporting)**
    *   目的：将上一步生成的空间表达特异性LD score与准备好的GWAS summary statistics进行关联分析，评估GWAS性状在特定空间区域/细胞类型中的富集程度，并生成报告。
    *   脚本：`05-spatial-ld&domain.sh` (生成任务列表), `submit_spatial_ldsc_array.sh` (提交任务)
    *   运行命令：
        ```bash
        # 生成任务列表 run_spatial_ldsc_jobs-hg.sh
        # 参数：工作目录，GWAS数据目录，并行数(提交脚本用)，内部线程数(gsmap用)，输出任务文件名，注释列名
        bash "05-spatial-ld&domain.sh" "/biostack/home/suchenyi/gsmap" "/biostack/home/suchenyi/graduate/gwas/gsmap_hg38" 10 50 "run_spatial_ldsc_jobs-hg.sh" "predicted_labels_plot" logs

        # 查看生成的任务数量
        wc -l run_spatial_ldsc_jobs-hg.sh

        # 提交任务 (示例：并行数为 1，日志输出到 logs/ 目录)
        # 注意：实际并行数取决于服务器资源和任务数量
        mkdir -p logs # 创建日志目录
        sbatch --array=1-1000%10 submit_spatial_ldsc_array.sh
        ```
    *   **监控与错误处理**：
        *   使用 `pgrep -f gsmap` 或 `htop` 监控 `gsmap spatial_ldsc` 或 `gsmap report` 进程。
        *   **检查已完成的任务**：(假设gsmap的日志输出在`logs/`目录下，且包含"Finished"和"report"字样)
          ```bash
          grep -H "Finished" nohup_logs/*.log | grep report | cut -d":" -f1 | sed 's|nohup_logs/||g' > success_hg-05.txt
          ```
        *   **生成重跑脚本**：如果任务因中断或错误未全部完成，可以比较任务列表和成功列表，生成未完成任务的脚本。
          ```bash
          grep -v -f success_hg-05.txt run_spatial_ldsc_jobs-hg.sh > 05-run_remain-hg.sh
          wc -l success_hg-05.txt 05-run_remain-hg.sh
          # 重新提交未完成的任务
          ```
        *   **强制停止**：如果需要停止所有 `gsmap` 进程，可以使用 `pkill -f gsmap` (谨慎使用)。