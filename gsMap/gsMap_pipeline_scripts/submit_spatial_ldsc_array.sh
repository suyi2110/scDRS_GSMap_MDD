#!/bin/bash
#SBATCH -J spatial_ldsc
#SBATCH -p cpu_batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=75G
#SBATCH -t 24:00:00
#SBATCH -o nohup_logs/%A_%a.out
#SBATCH -e nohup_logs/%A_%a.err

TASK_FILE="run_spatial_ldsc_jobs_human.sh"

CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $TASK_FILE)

echo "Running task ${SLURM_ARRAY_TASK_ID}"
echo "$CMD"

eval $CMD
