#!/bin/bash

#SBATCH --array=0
# #SBATCH --array=0-95
#SBATCH --time=1-0

#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/summarize_%A_%a"

export PYTHONPATH=/home/tcomi/projects/aclark4_introgression/code/

module load anaconda3
conda activate introgression3

#SLURM_ARRAY_TASK_ID=0
ARGS="_test .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"

python ${PYTHONPATH}analyze/summarize_region_quality_main.py $SLURM_ARRAY_TASK_ID $ARGS
