#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=2 # each task uses 1 cpu
#SBATCH --partition=urtgen_unlimit
#SBATCH --mem-per-cpu=140000 #200GB


source /home/u/f058977/miniconda3/etc/profile.d/conda.sh
conda activate MAGMA_classification_env
# alias R="R --no-save --no-restore --quiet"
R --no-save --no-restore CMD BATCH ./Script_Classification_v2/CALLING_SCRIPT_v2.R


