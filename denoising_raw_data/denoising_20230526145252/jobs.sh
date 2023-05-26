#!/bin/bash
  
#SBATCH --job-name=Array
#SBATCH --time=24:00:00
#SBATCH --qos=1day
#SBATCH --mem=1G
#SBATCH --output=/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/denoising_20230526145252/results_%a
#SBATCH --error=/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/denoising_20230526145252/errors_%a

source $HOME/.bashrc
conda activate 

JOBSFILES=/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/denoising_20230526145252/commands.cmd
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $JOBSFILES)
eval $SEED
