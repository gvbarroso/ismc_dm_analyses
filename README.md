# Quantifying the determinants of genome-wide distribution of diversity in the fruit fly
R scripts for reproducing results from Barroso & Dutheil 2020

## Tool scripts
simulate_seqs.R => simulates rho and theta landscapes and performs coalescent simulations with SCRM

bin_sim_maps.R => bins simulated rho and theta landscapes in different window sizes

bin_tmrca.R => bins simulated TMRCA landscapes in different window sizes

compute_pi.R => computes Tajima's pi for sequences simulated with simulate_seqs.R and bins them in different window sizes

## Data analyses scripts
R2_evol_sims.R => computes variance explained by each genomic landscape in the 12 evolutionary scenarios explored in our simulated studies

dm_analyses.Rmd => performs all analyses in Drosophila data and Drosophila-like simulations

To reproduce the results from the manuscript, data must be first downloaded from FigShare:

10.6084/m9.figshare.13164320

and extracted using the command
```bash
tar xvJf root.tar.gz
```

Once extracted, change into the "root" directory and run dm_analyses.Rmd to generate a notebook in PDF format with all intermediate and final results from Drosophila melanogaster analyses (both real and simulated data). Generating this notebook takes < 2 min. Run R2_evol_sims.R to generate R^2 plots for our simulation study of alternative evolutionary scenarios
