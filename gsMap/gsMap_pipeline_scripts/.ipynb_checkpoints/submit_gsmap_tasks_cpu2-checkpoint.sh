#!/bin/bash

MAX_JOBS=$2
INPUT_FILE="$1"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file $INPUT_FILE does not exist." >&2
    exit 1
fi

while IFS= read -r job; do
    while true; do
        numJobs=$(pgrep -cx gsmap)
        if [ "$numJobs" -lt "$MAX_JOBS" ]; then
            break
        fi
        sleep 3
    done

    echo "Submitting job: $job"
    eval "$job" &   # 关键：后台运行
    sleep 1         # 给系统一点缓冲
done < "$INPUT_FILE"

echo "All tasks in $INPUT_FILE have been submitted."

# 等待所有 gsmap 完成
while pgrep -x gsmap > /dev/null; do
    sleep 10
done

echo "All tasks in $INPUT_FILE have completed."
