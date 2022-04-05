#!/bin/bash
#SBATCH -p medium
#SBATCH -c 10
#SBATCH --mem-per-cpu=32gb 
#SBATCH -o scriptrmarkdown-japonica_2_xgbreg_phr
#SBATCH --time=48:00:00

module load r/4.0.3
source activate lightgbm_old

Rscript --max-ppsize=500000 script_forward_japonica_stacking_4_GY.R
