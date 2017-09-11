#!/bin/bash -l

#SBATCH
#SBATCH --job-name=ihiinc
#SBATCH --time=2:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=shared
#SBATCH --mail-type=end
#SBATCH

module load R

Rscript Xpert_Ultra_run_model_plosmed.R India NA highinc
# inputs to R script are setting, fixed params, tag ("overall" to run main model, or alternatives can specify sensitivity analyses, etc)
