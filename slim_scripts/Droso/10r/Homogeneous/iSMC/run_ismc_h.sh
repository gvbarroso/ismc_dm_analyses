#!/bin/bash
#SBATCH --job-name=dm2L_10r_h_ismc_h
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=400:00:00
#SBATCH --mem=100G
#SBATCH --error=dm2L_10r_h_ismc_h.%J.err
#SBATCH --output=dm2L_10r_h_ismc_h.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dutheil@evolbio.mpg.de
#SBATCH --partition=global

run_ismc() {
  mkdir -p Homogeneous/rep$1
  cd Homogeneous/rep$1
  $HOME/.local/bin/ismc params=../../ismc_h.bpp DATA=dm2L_bgs_rep${1}_5
  rm -rf ziphmm*
  cd ../..
}

run_ismc $SLURM_ARRAY_TASK_ID

echo "Done."

