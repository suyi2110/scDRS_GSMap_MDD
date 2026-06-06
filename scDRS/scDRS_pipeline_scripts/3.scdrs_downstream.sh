#!/bin/bash
#SBATCH -J scdrs_down
#SBATCH -p cpu_batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=80G
#SBATCH -t 72:00:00
#SBATCH --array=0-89%10
#SBATCH -o logs/scdrs_%A_%a.out
#SBATCH -e logs/scdrs_%A_%a.err

# 激活环境
source ~/.bashrc
module load Anaconda/
conda activate scdrs_env

# 读取 trait 名
TRAIT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" female_score/traits_list.txt)

echo "Running trait: $TRAIT"

scdrs perform-downstream \
    --h5ad-file /biostack/home/suchenyi/graduate/scdrs/Annotated_Female1.h5ad \
    --score-file female_score1/${TRAIT}.full_score.gz \
    --out-folder ./female_downstream1/ \
    --group-analysis "cell_type" \
    --flag-filter-data False \
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /biostack/home/suchenyi/graduate/scdrs/Annotated_Female1.h5ad \
    --score-file female_score1/${TRAIT}.full_score.gz \
    --out-folder ./female_downstream1/ \
    --group-analysis "cell_subtype" \
    --flag-filter-data False \
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /biostack/home/suchenyi/graduate/scdrs/Annotated_Female1.h5ad \
    --score-file female_score1/${TRAIT}.full_score.gz \
    --out-folder ./female_downstream1/ \
    --group-analysis "region_con" \
    --flag-filter-data False \
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /biostack/home/suchenyi/graduate/scdrs/Annotated_Female1.h5ad \
    --score-file female_score1/${TRAIT}.full_score.gz \
    --out-folder ./female_downstream1/ \
    --group-analysis "ct_con" \
    --flag-filter-data False \
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /biostack/home/suchenyi/graduate/scdrs/Annotated_Female1.h5ad \
    --score-file female_score1/${TRAIT}.full_score.gz \
    --out-folder ./female_downstream1/ \
    --group-analysis "cts_con" \
    --flag-filter-data False \
    --flag-raw-count True
