Run iSMC

To save time, we first estimate demographic parameters with a homogeneous model:
```bash
run_ismc() {
  mkdir -p Homogeneous/rep$1
  cd Homogeneous/rep$1
  $HOME/.local/bin/ismc params=../../ismc_h.bpp DATA=dm2L_bgs_rep${1}_5
  rm -rf ziphmm*
  cd ../..
}

run_ismc 1

```
Or on the cluster:
```
sbatch run_ismc_h.sh
```
