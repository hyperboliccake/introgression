#!/bin/bash

#SBATCH --time=0-4
#SBATCH --array=0,2-5
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/summarize_%A_%a"

export PYTHONPATH=/home/tcomi/projects/aclark4_introgression/code/

module load anaconda3
conda activate introgression3

ARGS="_test .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"

python ${PYTHONPATH}analyze/summarize_region_quality_main.py $SLURM_ARRAY_TASK_ID $ARGS
