#!/bin/bash
#SBATCH -J scdrs_down
#SBATCH -p cpu_batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=80G
#SBATCH -t 72:00:00
#SBATCH --array=0-18
#SBATCH -o logs/scdrs_%A_%a.out
#SBATCH -e logs/scdrs_%A_%a.err

# 激活环境
source ~/.bashrc
module load Anaconda/
conda activate scdrs_env

# 读取 trait 名
TRAIT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" score_hg38/retry_list.txt)

echo "Running trait: $TRAIT"

scdrs perform-downstream \
    --h5ad-file /biostack/home/suchenyi/scdrs/raw_hg.h5ad \
    --score-file score_hg38/${TRAIT}.full_score.gz \
    --out-folder ./downstream_hg38/ \
    --group-analysis "section" \
    --flag-filter-data False \
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /biostack/home/suchenyi/scdrs/raw_hg.h5ad \
    --score-file score_hg38/${TRAIT}.full_score.gz \
    --out-folder ./downstream_hg38/ \
    --group-analysis "condition" \
    --flag-filter-data False \
    --flag-raw-count True

