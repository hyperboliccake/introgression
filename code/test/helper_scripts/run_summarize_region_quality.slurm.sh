#!/bin/bash

#SBATCH --time=0-4

#SBATCH --array=0-5
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/summarize_%A_%a"

export PYTHONPATH=/home/tcomi/projects/aclark4_introgression/code/

module load anaconda3
conda activate introgression3

ARGS="p4e2 .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"

if [[ $SLURM_ARRAY_TASK_ID == 0 ]]; then
    start=10
else
    start=0
fi

start=$(($start + $SLURM_ARRAY_TASK_ID * 16))
end=$(($SLURM_ARRAY_TASK_ID * 16 +15))

for id in $(seq $start $end); do
    echo $id 
    python ${PYTHONPATH}analyze/summarize_region_quality_main.py $id $ARGS > /dev/null
done
