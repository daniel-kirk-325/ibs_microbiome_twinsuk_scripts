#!/bin/bash -l
#SBATCH --job-name=phyl
#SBATCH --output=log.out
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --array=0-15

module load python/3.11.6-gcc-13.2.0
source venv/bin/activate
python script.py
deactivate