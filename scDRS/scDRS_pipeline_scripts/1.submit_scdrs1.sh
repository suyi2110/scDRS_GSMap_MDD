#!/bin/bash
#SBATCH --job-name=scDRS_computescore          # 任务名
#SBATCH --partition=cpu_batch           
#SBATCH --nodes=1                   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G                     
#SBATCH --output=res_%j.log
#SBATCH --error=err_%j.log

# 1. 环境加载
source $(conda info --base)/etc/profile.d/conda.sh
conda activate scdrs_env

# 2. 路径设置
WORKING_DIR="/biostack/home/suchenyi/graduate/scdrs"
#女性scDRS输入数据
H5AD_FILE="/biostack/home/suchenyi/graduate/scdrs/Annotated_Female1.h5ad"
#男性scDRS输入数据
#H5AD_FILE="/biostack/home/suchenyi/scdrs/raw_hg_combined.h5ad"
#基因集文件
GS_FILE="/biostack/home/suchenyi/graduate/gwas/scdrs_hg38/geneset.gs"
OUT_DIR="${WORKING_DIR}/female_score1"

mkdir -p $OUT_DIR
cd $WORKING_DIR

# 3. 执行计算
scdrs compute-score \
    --h5ad-file $H5AD_FILE \
    --h5ad-species human \
    --gs-file $GS_FILE \
    --gs-species human \
    --flag-filter-data True \
    --flag-raw-count True \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --n-threads 16 \
    --out-folder $OUT_DIR

echo "Task Completed: `date`"
