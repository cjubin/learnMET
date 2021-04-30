#!/bin/bash
#SBATCH -p medium
#SBATCH -c 10
#SBATCH --mem-per-cpu=18gb 
#SBATCH -o scriptrmarkdown-%J
#SBATCH --array=2
#SBATCH --time=48:00:00


Rscript --max-ppsize=500000 example_workflow.R