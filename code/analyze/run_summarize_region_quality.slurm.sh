#!/bin/bash
#SBATCH --array=0-95
#SBATCH --time=1-0 --qos=10hr
#SBATCH -n 1

ARGS=$(head -n 1 predict_args.txt)

python summarize_region_quality_main.py $SLURM_ARRAY_TASK_ID $ARGS
