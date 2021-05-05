#!/bin/bash
#SBATCH -c 8  # Number of Cores per Task
#SBATCH -p cee_water_cjgleason
#SBATCH --mem=25000  # Requested Memory
#SBATCH -t 500:00:00  # Job time limit
#SBATCH -o out.txt  # %j = job ID
#SBATCH -e err.txt

source ~/.bashrc

module load R
module load gcc/9.3.0

Rscript src/main.R
