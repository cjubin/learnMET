#!/bin/bash
#SBATCH -p fat
#SBATCH -c 3
#SBATCH -N 1
#SBATCH --mem-per-cpu=25gb 
#SBATCH -o scriptrmarkdown-%J
#SBATCH --array=1-4
#SBATCH -t 0-32:00:00

module load anaconda3/2021.05
source activate lightgbm_old
Rscript --max-ppsize=500000 -e "rmarkdown::render('BGLR_G_E_W_GW_GY.Rmd','html_document',params=list(sets_predictors = 'G+WC+SC+Lon+Lat',env_to_predict=(as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')))),output_file=paste0(as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')),'default_parameters','_G+WC+SC+Lon+Lat','.html'))"


