Run iSMC

To save time, we first estimate demographic parameters with a homogeneous model:
```bash
sbatch run_ismc_h.sh
```

We then fit the heterogeneous model:
```bash
sbatch run_ismc.sh
```

And get the posterior decoding:
```bash
sbatch run_ismc_decoding.sh
```

Then compute averages in windows:
```bash
for REP in {1..10}; do
  echo $REP
  cd NonHomogeneous/rep$REP
  $HOME/.local/bin/ismc_mapper params=../../ismc_mapper.bpp DATA=dm2L_bgs_rep${REP}_5 bin_rate=rho fasta_seqs=true tmrca=true
  $HOME/.local/bin/ismc_mapper params=../../ismc_mapper.bpp DATA=dm2L_bgs_rep${REP}_5 bin_rate=theta fasta_seqs=false tmrca=false
  cd ../..
done
```
