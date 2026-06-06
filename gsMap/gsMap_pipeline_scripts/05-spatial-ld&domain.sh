#!/bin/bash
# ===========================================================
# gsmap run_spatial_ldsc + run_report 批量命令生成脚本（使用样本目录下 _annotated.h5ad）
# 用法: bash generate_spatial_ldsc_jobs.sh <WORKDIR> <GWAS_DIR> <NUM_PROCESSES> <TOP_CORR_GENES> <OUTPUT_SCRIPT> <ANNOT> <LOG_DIR>
# ===========================================================

WORKDIR=$1              # 工作目录，例如 "work-prior"
GWAS_DIR=$2             # GWAS 性状文件目录，例如 "data/GWAS/Traits-1"
NUM_PROCESSES=$3        # 并行进程数，例如 22
TOP_CORR_GENES=$4       # 报告中显示的相关基因数，例如 30
OUTPUT_SCRIPT=$5        # 输出脚本文件名，例如 "run_spatial_ldsc_jobs.sh"
ANNOT=$6                # 注释列，例如 "predicted_labels_plot"
LOG_DIR=$7              # 日志目录，例如 "nohup_logs"

# LDSC 权重路径
W_FILE="/biostack/home/suchenyi/gsmap/gsMap_resource/LDSC_resource/weights_hm3_no_hla/weights."

# 创建日志目录
mkdir -p "$LOG_DIR"
> "$OUTPUT_SCRIPT"  # 清空输出脚本文件

# 遍历 WORKDIR 下每个样本目录
for sample_dir in "$WORKDIR"/*_rep*; do
    sample_name=$(basename "$sample_dir")
    h5ad_file="$sample_dir/${sample_name}_annotated.h5ad"

    # 检查 h5ad 是否存在
    if [[ ! -f "$h5ad_file" ]]; then
        echo "[跳过] h5ad 文件不存在: $h5ad_file"
        continue
    fi

    # 遍历每个 GWAS 性状文件
    for trait_file in "$GWAS_DIR"/*.sumstats.gz; do
        trait_name=$(basename "$trait_file" .sumstats.gz | sed 's/trait_//')

        # spatial_ldsc 命令
        cmd_spatial_ldsc="gsmap run_spatial_ldsc \
            --workdir '$WORKDIR' \
            --sample_name '$sample_name' \
            --trait_name '$trait_name' \
            --sumstats_file '$trait_file' \
            --w_file '$W_FILE' \
            --num_processes $NUM_PROCESSES"

        # report 命令
        cmd_report="gsmap run_report \
            --workdir '$WORKDIR' \
            --sample_name '$sample_name' \
            --trait_name '$trait_name' \
            --annotation '$ANNOT' \
            --sumstats_file '$trait_file' \
            --top_corr_genes $TOP_CORR_GENES"

        # 合并命令（不加 &，提交脚本控制并发）
        combined_cmd="$cmd_spatial_ldsc"

        # 构造日志
        log_file="$LOG_DIR/${sample_name}_${trait_name}.log"

        # 写入输出脚本
        echo "nohup bash -c \"$combined_cmd\" > \"$log_file\" 2>&1" >> "$OUTPUT_SCRIPT"
    done
done

echo "所有命令已生成完毕，保存在 $OUTPUT_SCRIPT，总计 $(wc -l < $OUTPUT_SCRIPT) 行。"
