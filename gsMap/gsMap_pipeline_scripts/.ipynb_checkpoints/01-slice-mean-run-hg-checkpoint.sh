#!/bin/bash

# 设置数据目录和输出目录
DATA_DIR=$1  # "data/prior-stereoseq"
OUTPUT_DIR=$2 # "./work-prior"
# HOMOLOG_FILE=$3 # "data/mouse_human_homologs.txt"
LOG_DIR="nohup_logs"

# 确保输出目录和日志目录存在
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# 获取所有_rep1和_rep2的文件
files=$(ls "$DATA_DIR"/*rep[12]_*)

# 使用关联数组存储brain_region_condition组
declare -A groups

# 遍历文件并分组
for file in $files; do
    # 提取brain_region_condition，例如DRN_Control
    base=$(basename "$file" .h5ad)
    brain_region_condition=$(echo "$base" | sed 's/_rep[12]_reannotated//')
    # 将文件添加到对应的组
    groups["$brain_region_condition"]+="$file "
done

# 遍历每个brain_region_condition组并生成命令
for brc in "${!groups[@]}"; do
    # 获取该组的文件列表
    file_list=(${groups["$brc"]})
    # 确保组中有两个文件
    if [ ${#file_list[@]} -eq 2 ]; then
        # 提取样本名称
        sample1=$(basename "${file_list[0]}" _reannotated.h5ad)
        sample2=$(basename "${file_list[1]}" _reannotated.h5ad)
        # 生成输出文件名
        output_file="$OUTPUT_DIR/${brc}_slice_mean.parquet"
        # 生成gsmap命令
        cmd="gsmap create_slice_mean \
            --sample_name_list '$sample1' '$sample2' \
            --h5ad_list '${file_list[0]}' '${file_list[1]}' \
            --slice_mean_output_file '$output_file' \
            --data_layer 'count'"

        # 使用nohup后台运行命令，并将日志保存到LOG_DIR
        echo "正在提交任务: $brc"
        nohup bash -c "$cmd" > "$LOG_DIR/${brc}.log" 2>&1 &

        # 控制并发任务数量，保持10个任务同时运行
        while true; do
            numJobs=$(pgrep -cx gsmap)
            if [ "$numJobs" -lt 30 ]; then
                break
            fi
            sleep 3
        done
    else
        echo "警告: $brc 组的文件数量不为2: ${file_list[@]}"
    fi
done

# 脚本结束提示
echo "所有任务已提交完毕。"



