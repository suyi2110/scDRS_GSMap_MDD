#!/bin/bash
# =======================
# gsmap run_find_latent_representations 脚本（修改版）
# =======================

# 设置相关目录和参数
WORKDIR="/biostack/home/suchenyi/graduate/gsmap/human_h5ad"  #"work-prior"                      # 工作目录
INPUT_H5AD_DIR="/biostack/home/suchenyi/graduate/gsmap/data"   #"data/prior-stereoseq"     # 输入 h5ad 文件目录
ANNOTATION="manual_cell_type_reannotated"   #"class_level"                   # 注释字段名称
DATA_LAYER="count"                         # 数据层名称
OUTPUT_SCRIPT="/biostack/home/suchenyi/graduate/gsmap/run_find_latent_jobs-hg.sh"

# 确保日志目录存在
mkdir -p "nohup_logs"
rm -f "$OUTPUT_SCRIPT"

# 遍历输入目录下所有以 _gsmap.h5ad 结尾的文件
for file in "$INPUT_H5AD_DIR"/*_reannotated.h5ad; do
    # 从文件名中提取样本名（保留 _rep1 信息）
    sample=$(basename "$file" _reannotated.h5ad)  # 示例：DRN_Control_rep1

    # 构造 gsmap 命令
    cmd="gsmap run_find_latent_representations \
        --workdir '$WORKDIR' \
        --sample_name '$sample' \
        --input_hdf5_path '$file' \
        --data_layer 'count'\
        --annotation '$ANNOTATION'"

    echo "提交任务: $sample"
    # echo "$cmd"

    # 后台运行任务，日志输出至 nohup_logs 目录
    # nohup bash -c "$cmd" > "nohup_logs/${sample}.log" 2>&1 &
    # 输出将要执行的 nohup 命令（调试用）
    echo "nohup bash -c \"$cmd\" > \"nohup_logs/02-${sample}.log\" 2>&1 &"  >> "$OUTPUT_SCRIPT"
done

echo "所有任务已提交完毕。"
