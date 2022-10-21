#!/bin/bash
#SBATCH --job-name=dm2L_10r_nh_ismc_decoding
#SBATCH --array=1-2
#SBATCH --nodes=1
#SBATCH --ntasks=45
#SBATCH --time=1000:00:00
#SBATCH --mem=1500G
#SBATCH --error=dm2L_10r_nh_ismc_decoding.%J.err
#SBATCH --output=dm2L_10r_nh_ismc_decoding.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dutheil@evolbio.mpg.de
#SBATCH --partition=highmem

run_ismc() {
  cd NonHomogeneous/rep$1
  $HOME/.local/bin/ismc params=../../ismc_decoding.bpp DATA=dm2L_bgs_rep${1}_5
  rm -rf ziphmm*
  cd ../..
}

run_ismc $SLURM_ARRAY_TASK_ID

cd ..
echo "Done."

