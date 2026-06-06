#!/bin/bash
#SBATCH -J scdrs_plot
#SBATCH -p cpu_batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=20G
#SBATCH -t 24:00:00
#SBATCH -o logs/plot_%j.out
#SBATCH -e logs/plot_%j.err

source ~/.bashrc
module load Anaconda/
conda activate scdrs_env

python umap_sc.py
