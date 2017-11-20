#!/bin/bash

# By default it will run in the directory where sbatch was started.
# Or you can do something like this:
# cd $GRID_HOME

#SBATCH --array=1-17

# We need 20 GB for each instance (must be a good estimate).
# Note that SLURM assumes one core per task. For tasks using multiple
# threads, you should adjust --cpus-per-task and --mem-per-cpu.
# #SBATCH --mem-per-cpu=2048 --cpus-per-task=1
# Inform the scheduler we will finish within 2 minutes (optional).
# Formats: mins, hrs:mins, days-hrs:mins
# And use the 1hr QOS.
# #SBATCH --time=2 --qos=1hr
#SBATCH --time=0-20 --qos=1hr

# This specifies 32 tasks, which will require at least 2 nodes.
#SBATCH -n 1

# With 32 tasks allocated, srun will start 32 instances of
# the same program (in this case sh with an explicit command).
# You can use SLURM_PROCID to distinguish different tasks.
# It will have a value between 0 and 31 in this example.

ARGS=$(head -n $SLURM_ARRAY_TASK_ID predict_args.txt | tail -n 1)

python sim_predict_main.py sim_args.txt $ARGS
