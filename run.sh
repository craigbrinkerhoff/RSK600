#!/bin/bash
#SBATCH --job-name='BIKER'
#SBATCH -p cpu
#SBATCH -c 30  # Number of Cores per Task
#SBATCH --mem=64000  # Requested Memory
#SBATCH -t 5:00:00  # Job time limit
#SBATCH -o out.txt  # %j = job ID
#SBATCH -e err.txt

source ~/.bashrc

Rscript src/main.R
