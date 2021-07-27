#!/bin/bash
#SBATCH -p medium
#SBATCH -c 10
#SBATCH --mem-per-cpu=10gb 
#SBATCH -o scriptrmarkdown-analysis-indica-xgb
#SBATCH --time=48:00:00

module load r/4.0.3
#source activate lightgbm_old

Rscript --max-ppsize=500000 example_workflow.R