#!/bin/bash
#SBATCH -p medium
#SBATCH -c 10
#SBATCH --mem-per-cpu=18gb 
#SBATCH -o scriptrmarkdown-%J
#SBATCH --time=48:00:00

#module load anaconda3/2020.11 

#source activate lightgbm_old

Rscript --max-ppsize=500000 example_workflow.R