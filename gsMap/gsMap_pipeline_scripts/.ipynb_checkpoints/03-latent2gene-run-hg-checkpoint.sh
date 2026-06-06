#!/bin/bash
# =======================
# gsmap run_latent_to_gene 脚本（修改版）
# =======================

# 设置工作目录、数据目录、同源文件路径及其他参数
WORKDIR=$1          # "work-mm"                # gsmap 的工作目录，存放 slice_mean 文件的目录
INPUT_H5AD_DIR=$2   # "data/stereoseq/mm"      # 输入 h5ad 文件目录
# HOMOLOG_FILE=$3     # "data/mouse_human_homologs.txt"  # 同源基因文件
ANNOTATION=$3       # "class_level"            # 注释字段名称
NUM_NEIGHBOUR=51                    # 邻居数量
NUM_NEIGHBOUR_SPATIAL=201           # 空间邻居数量

# 确保日志目录存在
mkdir -p "nohup_logs"
rm "$4"

# 遍历输入目录下所有以 _gsmap.h5ad 结尾的文件
for file in "$INPUT_H5AD_DIR"/*_reannotated.h5ad; do
    # 从文件名中提取样本名（保留 _repX 信息）
    sample=$(basename "$file" _reannotated.h5ad)  # 示例：DRN_3D_rep1

    # 提取样本基础名称（去掉 _repX 后缀）
    base_sample=$(echo "$sample" | sed 's/_rep[12]\+$//')  # 示例：DRN_3D

    # 构造 slice_mean 文件路径
    slice_mean_file="$WORKDIR/${base_sample}_slice_mean.parquet"

    # 构造 gsmap run_latent_to_gene 命令
    cmd="gsmap run_latent_to_gene \
        --workdir '$WORKDIR' \
        --sample_name '$sample' \
        --annotation '$ANNOTATION' \
        --latent_representation 'latent_GVAE' \
        --num_neighbour $NUM_NEIGHBOUR \
        --num_neighbour_spatial $NUM_NEIGHBOUR_SPATIAL \
        --gM_slices '$slice_mean_file'"

    echo "提交任务: $sample"
    # echo "$cmd"

    # 输出 nohup 命令到文件，手动运行时使用
    echo "nohup bash -c \"$cmd\" > \"nohup_logs/03-${sample}.log\" 2>&1 &" >> "$4"
done

echo "所有任务已提交完毕。"
