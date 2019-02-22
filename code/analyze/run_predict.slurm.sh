#!/bin/bash

## SBATCH --array=1
#SBATCH --time=6-0

#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/predict_%A"

# ARGS=$(head -n $SLURM_ARRAY_TASK_ID predict_args.txt | tail -n 1)

export PYTHONPATH=/home/tcomi/projects/aclark4_introgression/code/

ARGS=$(head -n 1 predict_args.txt | tail -n 1)

module load anaconda

conda activate introgression

python predict_main.py $ARGS
