
Simulations
===========

Run 10 replicates on the cluster:
```bash
sbatch runSLiM.sh
```
This will stop after 700,000 generations.

Recapitate and sample tree sequence:
```bash
for i in {1..10}; do
  echo "Rep $i..."
  python3.9 recapitateAndSample.py dm2L_bgs_rep$i $((42 + i)) > dm2L_bgs_rep${i}_msprime_recapitation.log
done
```

Add mutations to the tree sequence. First with a uniform mutation rate:
```bash
mkdir Homogeneous
for i in {1..10}; do
  echo "Rep $i..."
  python3.9 ../simBGS.py dm2L_bgs_rep$i $((42 + i)) > Homogeneous/dm2L_bgs_rep${i}_msprime_mutation.log
done
```
Compress VCF files:
```bash
for i in {1..10}; do
  echo "Rep $i..."
  bgzip Homogeneous/dm2L_bgs_rep${i}_5.vcf
done
```

Then with a variable mutation rate.
```bash
mkdir NonHomogeneous
for i in {1..10}; do
  echo "Rep $i..."
  python3.9 ../simBGS_varmut.py dm2L_bgs_rep$i $((42 + i)) > NonHomogeneous/dm2L_bgs_rep${i}_msprime_mutation.log
done
```
Compress VCF files:
```bash
for i in {1..10}; do
  echo "Rep $i..."
  bgzip NonHomogeneous/dm2L_bgs_rep${i}_5.vcf
done
```

Compare with a purely neutral scenario (we still use SLiM so that all other things are equal):
```bash
sbatch runSLiMneutral.sh
```
Recapitate and sample tree sequence:
```bash
for i in {1..10}; do
  echo "Rep $i..."
  python3.9 recapitateAndSample.py dm2L_neutral_rep$i $((42 + i)) > dm2L_neutral_rep${i}_msprime_recapitation.log
done
```

Compute statistics
==================

Compute mean TMRCA in windows:
```bash
for i in {1..10}; do
  python3.9 ../getTmrca.py dm2L_bgs_rep$i $((42 + i)) 50000
  python3.9 ../getTmrca.py dm2L_bgs_rep$i $((42 + i)) 200000
  python3.9 ../getTmrca.py dm2L_bgs_rep$i $((42 + i)) 1000000
done
```

Compute diversity:
```bash
cd Homogeneous
for i in {1..10}; do
  vcftools --gzvcf dm2L_bgs_rep${i}_5.vcf.gz \
           --window-pi 50000 \
           --out dm2L_bgs_rep${i}_5_w50000
  vcftools --gzvcf dm2L_bgs_rep${i}_5.vcf.gz \
           --window-pi 200000 \
           --out dm2L_bgs_rep${i}_5_w200000
  vcftools --gzvcf dm2L_bgs_rep${i}_5.vcf.gz \
           --window-pi 1000000 \
           --out dm2L_bgs_rep${i}_5_w1000000
done
cd ..
cd NonHomogeneous
for i in {1..10}; do
  vcftools --gzvcf dm2L_bgs_rep${i}_5.vcf.gz \
           --window-pi 50000 \
           --out dm2L_bgs_rep${i}_5_w50000
  vcftools --gzvcf dm2L_bgs_rep${i}_5.vcf.gz \
           --window-pi 200000 \
           --out dm2L_bgs_rep${i}_5_w200000
  vcftools --gzvcf dm2L_bgs_rep${i}_5.vcf.gz \
           --window-pi 1000000 \
           --out dm2L_bgs_rep${i}_5_w1000000
done
cd ..
```

Merge all informations into one big data file:
```bash
mkdir RealLandscapes
```

```r
merge.data <- function(win, repl, output) {
  dat.rec <- read.table(paste0("../RecombinationMap_", win, "kb.csv"), sep = ",", header = TRUE)
  dat.rec$RecRate <- dat.rec$Rate * 10 #In these simulations, the original recombination rate was scaled by 10 to maintain the rho/theta ratio)
  dat.rec$Rate <- NULL
  dat.mut <- read.table(paste0("../MutationMap_", win, "kb.csv"), sep = ",", header = TRUE)
  dat.mut$MutRate <- dat.mut$Rate
  dat.mut$Rate <- NULL
  dat.tmrca <- read.table(paste0("dm2L_bgs_rep", repl, "_w", win, "000.csv"), sep = ",",  header = TRUE)
  dat.tmrca$Start <- dat.tmrca$Start + 1
  dat.tmrca$Stop <- dat.tmrca$Stop
  dat.h.div <- read.table(paste0("Homogeneous/dm2L_bgs_rep", repl, "_5_w", win, "000.windowed.pi"), sep = "\t", header = TRUE)
  dat.h.div$CHROM <- NULL
  names(dat.h.div) <- c("Start", "Stop", "NbVariants.Homogeneous", "Pi.Homogeneous")
  dat.nh.div <- read.table(paste0("NonHomogeneous/dm2L_bgs_rep", repl, "_5_w", win, "000.windowed.pi"), sep = "\t", header = TRUE)
  dat.nh.div$CHROM <- NULL
  names(dat.nh.div) <- c("Start", "Stop", "NbVariants.NonHomogeneous", "Pi.NonHomogeneous")
  dat <- merge(dat.rec, dat.mut, by = c("Start", "Stop"), sort = FALSE)
  dat <- merge(dat, dat.tmrca, by = c("Start", "Stop"), sort = FALSE)
  dat <- merge(dat, dat.h.div, by = c("Start", "Stop"), sort = FALSE)
  dat <- merge(dat, dat.nh.div, by = c("Start", "Stop"), sort = FALSE)
  write.table(dat, file = output, sep = ",", row.names = FALSE, quote = FALSE)
}

for (sim in c("50", "200", "1000")) {
  for (repl in 1:10) {
    merge.data(sim, repl, paste0("RealLandscapes/RealLandscapes_", sim, "kb_rep", repl, ".csv"))
  }
}
```

Get TMRCA in the neutral case for comparison:
```bash
for i in {1..10}; do
  python3.9 ../getTmrca.py dm2L_neutral_rep$i $((42 + i)) 50000
  python3.9 ../getTmrca.py dm2L_neutral_rep$i $((42 + i)) 200000
  python3.9 ../getTmrca.py dm2L_neutral_rep$i $((42 + i)) 1000000
done
```

