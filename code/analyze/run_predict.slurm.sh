#!/bin/bash
#SBATCH --array=1
#SBATCH --time=6-0 --qos=10hr
#SBATCH -n 1

ARGS=$(head -n $SLURM_ARRAY_TASK_ID predict_args.txt | tail -n 1)

python predict_main.py $ARGS
