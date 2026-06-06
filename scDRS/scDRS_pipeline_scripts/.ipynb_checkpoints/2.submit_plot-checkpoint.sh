#!/bin/bash
#SBATCH --job-name=scDRS_plot
#SBATCH --partition=cpu_batch
#SBATCH --nodes=2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --mem=10G 
#SBATCH --time=05:00:00      
#SBATCH --output=plot_%j.log
#SBATCH --error=plot_%j.err

# 环境激活
source $(conda info --base)/etc/profile.d/conda.sh
conda activate scdrs_env

# 执行脚本
python umap_plot.py
