#!/bin/bash
#SBATCH --job-name=scDRS_computescore          # 任务名
#SBATCH --partition=cpu_batch             # 请根据实际情况填写，如有GPU需求可填gpu
#SBATCH --nodes=1                       # 单节点运行即可，scDRS主要靠单机多线程
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G                      # 512G总内存，分给这个任务200G非常稳妥
#SBATCH --output=res_%j.log
#SBATCH --error=err_%j.log

# 1. 环境加载
source $(conda info --base)/etc/profile.d/conda.sh
conda activate scdrs_env

# 2. 路径设置 (建议全部使用绝对路径)
WORKING_DIR="/biostack/home/suchenyi/graduate/scdrs"
H5AD_FILE="/biostack/home/suchenyi/scdrs/raw_hg.h5ad"
GS_FILE="/biostack/home/suchenyi/graduate/gwas/scdrs/geneset3.gs"
OUT_DIR="${WORKING_DIR}/trait_analysis3"

mkdir -p $OUT_DIR
cd $WORKING_DIR

# 3. 执行计算
# 增加 --n-threads 参数以匹配 SBATCH 分配的 CPU 核心数
scdrs compute-score \
    --h5ad-file $H5AD_FILE \
    --h5ad-species human \
    --gs-file $GS_FILE \
    --gs-species human \
    --flag-filter-data True \
    --flag-raw-count True \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --n-threads 48 \
    --out-folder $OUT_DIR

echo "Task Completed: `date`"
