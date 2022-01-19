#!/bin/bash

#SBATCH --job-name=bgs-slim10
#SBATCH --ntasks=1 # num of cores
#SBATCH --nodes=1
#SBATCH --time=300:00:00 # in hours
#SBATCH --mem=100G 
#SBATCH --error=log/bgs_slim10.%J.err
#SBATCH --output=log/bgs_slim10.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=####
#SBATCH --partition=global 
#SBATCH --array=1-10

~/.local/bin/slim -d seed=$((SLURM_ARRAY_TASK_ID + 42)) -d rep=$SLURM_ARRAY_TASK_ID Simulate.slim
  
