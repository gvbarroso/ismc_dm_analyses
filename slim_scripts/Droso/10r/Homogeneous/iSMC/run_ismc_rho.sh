#!/bin/bash
#SBATCH --job-name=dm2L_10r_nh_ismc_rho
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=1000:00:00
#SBATCH --mem=256G
#SBATCH --error=dm2L_10r_nh_ismc_rho.%J.err
#SBATCH --output=dm2L_10r_nh_ismc_rho.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dutheil@evolbio.mpg.de
#SBATCH --partition=standard

run_ismc() {
  mkdir -p NonHomogeneousRho/rep$1
  cp Homogeneous/rep$1/dm2L_bgs_rep${1}_5_backup_params.txt NonHomogeneousRho/rep$1/
  cd NonHomogeneousRho/rep$1
  $HOME/.local/bin/ismc params=../../ismc_rho.bpp DATA=dm2L_bgs_rep${1}_5 resum_optim=true
  rm -rf ziphmm*
  cd ../..
}

run_ismc $SLURM_ARRAY_TASK_ID

echo "Done."

