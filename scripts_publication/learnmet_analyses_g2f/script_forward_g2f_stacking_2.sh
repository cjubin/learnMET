#!/bin/bash
#SBATCH -p fat
#SBATCH -c 10
#SBATCH --mem-per-cpu=32gb 
#SBATCH -o scriptrmarkdown-g2f_2_stacking
#SBATCH --time=48:00:00

module load r/4.0.3
source activate lightgbm_old

Rscript --max-ppsize=500000 script_forward_g2f_stacking_2.R
