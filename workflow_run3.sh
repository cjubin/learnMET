#!/bin/bash
#SBATCH -p medium
#SBATCH -c 4
#SBATCH --mem-per-cpu=10gb 
#SBATCH -o scriptrmarkdown-analysis-indica-svm
#SBATCH --time=48:00:00

#module load anaconda3/2020.11 

#source activate lightgbm_old

Rscript --max-ppsize=500000 example_workflow3.R