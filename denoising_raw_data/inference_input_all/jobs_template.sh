#!/bin/bash
  
#SBATCH --job-name=Array
#SBATCH --time=24:00:00
#SBATCH --qos=1day
#SBATCH --mem=1G
#SBATCH --output=path_to_output
#SBATCH --error=path_to_error

source $HOME/.bashrc
conda activate 

JOBSFILES=path_to_commands
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $JOBSFILES)
eval $SEED
