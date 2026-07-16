#!/bin/bash -l
#SBATCH --job-name=NA
#SBATCH --output=logs_perm/log_perm_%a.out
#SBATCH --time=07:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --array=1-999


module load anaconda3/2022.10-gcc-13.2.0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /cephfs/volumes/hpc_data_usr/k2365231/7dce5135-5ef3-4ed3-9176-23b8192716f8/2026/network_analysis/r_env

export I=${SLURM_ARRAY_TASK_ID}

echo "R version:"
R --version

Rscript network_permutations.R