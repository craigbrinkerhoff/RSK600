#!/bin/bash
#SBATCH --job-name='BIKER'
#SBATCH -p cpu
#SBATCH -c 48  # Number of Cores per Task
#SBATCH --mem=128000  # Requested Memory
#SBATCH -t 5:00:00  # Job time limit
#SBATCH -o out.txt  # %j = job ID
#SBATCH -e err.txt

source ~/.bashrc

module load R
module load gcc/9.3.0

Rscript src/main.R
