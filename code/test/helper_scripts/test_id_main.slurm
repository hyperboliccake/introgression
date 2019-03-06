#!/bin/bash

## SBATCH --array=1
#SBATCH --time=2-0

#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/id_%A"

# ARGS=$(head -n $SLURM_ARRAY_TASK_ID predict_args.txt | tail -n 1)

export PYTHONPATH=/home/tcomi/projects/aclark4_introgression/code/

#ARGS="p4e2 .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"
ARGS="_test .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"

#Make sure chrms is set to only I
#ARGS="_chr1_test .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"
#ARGS="_chr1 .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"

module load anaconda3
conda activate introgression3

python $PYTHONPATH/analyze/id_regions_main.py $ARGS
