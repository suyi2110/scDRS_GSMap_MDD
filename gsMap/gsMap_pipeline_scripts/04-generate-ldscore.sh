#!/bin/bash
# =======================
# gsmap run_generate_ldscore 脚本（修改版，每个任务使用10个核心）
# =======================

# 设置相关目录和参数
WORKDIR=$1              # 工作目录，例如 "work-prior"
INPUT_H5AD_DIR=$2       # 输入 h5ad 文件目录，例如 "data/stereoseq/mm"
GTF_FILE=$3             # 例如 "data/hg_genome/gencode.v31lift37.annotation.gtf"
BFILE_ROOT="/biostack/home/suchenyi/graduate/gwas/hg38/reference/LDreference/chrom/1000G.EUR.QC"
KEEP_SNP_ROOT="/biostack/home/suchenyi/gsmap/gsMap_resource/LDSC_resource/hapmap3_snps/hm"
GENE_WINDOW_SIZE=50000
CHUNK_SIZE=500         # CELL
LOG_DIR="nohup_logs"   # 日志目录

# 确保日志目录存在
mkdir -p "$LOG_DIR"
rm -f "$4"  # 清空任务命令输出文件

# 获取总的核心数
total_cores=$(nproc)
block_size=10  # 每个任务使用10个核心
# 初始化任务计数器
task_index=0

# 遍历输入目录下所有以 _gsmap.h5ad 结尾的文件
for file in "$INPUT_H5AD_DIR"/*_reannotated.h5ad; do
    # 从文件名中提取样本名（保留 _repX 信息）
    sample=$(basename "$file" _reannotated.h5ad)  # 示例：DRN_3D_rep1

    # 遍历 1~22 号染色体
    for CHROM in {1..22}; do
        # 构造 gsmap 命令
        cmd="gsmap run_generate_ldscore \
            --workdir '$WORKDIR' \
            --sample_name '$sample' \
            --chrom $CHROM \
            --bfile_root '$BFILE_ROOT' \
            --keep_snp_root '$KEEP_SNP_ROOT' \
            --gtf_annotation_file '$GTF_FILE' \
            --gene_window_size '$GENE_WINDOW_SIZE' \
            --spots_per_chunk '$CHUNK_SIZE'"

        echo "提交任务: $sample 染色体 $CHROM"

        # 计算本任务分配的核心块（10个核心）
        start=$(( (task_index * block_size) % total_cores ))
        end=$(( start + block_size - 1 ))
        if [ $end -ge $total_cores ]; then
            # 如果超出最高核心，则从头开始补充
            part1=$(seq $start $((total_cores - 1)) | paste -sd, -)
            part2=$(seq 0 $((end - total_cores)) | paste -sd, -)
            cores="$part1,$part2"
        else
            cores=$(seq $start $end | paste -sd, -)
        fi

        # nohup 后台运行任务，使用 taskset 指定这10个核心
        # nohup taskset -c "$cores" bash -c "$cmd" > "$LOG_DIR/04-${sample}_chr${CHROM}.log" 2>&1 &
        # 将执行的命令写入任务日志文件（$4 指定的文件）
        # echo "nohup taskset -c \"$cores\" bash -c \"$cmd\" > \"${LOG_DIR}/04-${sample}_chr${CHROM}.log\" 2>&1 &" >> "$4"
        # echo "nohup bash -c \"$cmd\" > \"${LOG_DIR}/04-${sample}_chr${CHROM}.log\" 2>&1 &" >> "$4"
        # 0-45, 用于每次只跑一次
        taskset -c "$cores"
        echo "nohup taskset -c \"$cores\" bash -c \"$cmd\" > \"${LOG_DIR}/04-${sample}_chr${CHROM}.log\" 2>&1 &" >> "$4"

        # 任务计数器自增
        task_index=$((task_index + 1))
    done
done

echo "所有任务已提交完毕。"
