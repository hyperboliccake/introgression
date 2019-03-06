#!/bin/bash

## SBATCH --array=1
#SBATCH --time=6-0

#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/predict_%A"

# ARGS=$(head -n $SLURM_ARRAY_TASK_ID predict_args.txt | tail -n 1)

export PYTHONPATH=/home/tcomi/projects/aclark4_introgression/code/

#Make sure chrms is set to only I
#ARGS="_chr1_test .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"
ARGS="_test .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"

module load anaconda3
conda activate introgression3

#python $PYTHONPATH/analyze/predict_main.py $ARGS
# gzip after
gzip /tigress/tcomi/aclark4_temp/results/analysis_test/positions__test.txt
gzip /tigress/tcomi/aclark4_temp/results/analysis_test/probs__test.txt
