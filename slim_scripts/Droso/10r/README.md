
Simulations
===========

Run 10 replicates on the cluster:
```bash
sbatch runSLiM.sh
```
This will stop after 700,000 generations.
 
Add mutations to the tree sequence. First with a uniform mutation rate:
```bash
mkdir Homogeneous
for i in {1..10}; do
  echo "Rep $i..."
  python3.9 ../simBGS.py dm2L_bgs_rep$i $((42 + i)) > Homogeneous/dm2L_bgs_rep${i}_msprime.log
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
  python3.9 ../simBGS_varmut.py dm2L_bgs_rep$i $((42 + i)) > NonHomogeneous/dm2L_bgs_rep${i}_msprime.log
done
```
Compress VCF files:
```bash
for i in {1..10}; do
  echo "Rep $i..."
  bgzip NonHomogeneous/dm2L_bgs_rep${i}_5.vcf
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


Perform analyses
================

```r
dat <- read.table("dm2L_bgs_rep1_5_w50000.windowed.pi", header = TRUE)
taj <- read.table("dm2L_bgs_rep1_5_w50000.Tajima.D", header = TRUE)
snp <- read.table("dm2L_bgs_rep1_5_w50000.snpden", header = TRUE)
sin <- read.table("dm2L_bgs_rep1_5_w50000_maxmac1.snpden", header = TRUE)
m30 <- read.table("dm2L_bgs_rep1_5_w50000_maf30.snpden", header = TRUE)
rec <- read.table("Comeron_tables/Comeron_100kb_chr2L.txt", header = F)
tmrca <- read.table("dm2L_bgs_rep1_w50000.tmrca_sample5.csv", header = T, sep = ",")
exons <- read.table("GeneAnnotations/exons_chr2L.csv", sep = ",")
names(exons) <- c("xmin", "xmax")

# Compute gene density:
library(IRanges)
ranges.exons <- IRanges(start = exons$xmin, end = exons$xmax)
s <- seq(1, 23e6, by = 50000)
dens <- numeric(length(s))
for (i in 1:length(s)) {
  r <- IRanges(start = s[i], end = s[i] + 50e3)
  q <- intersect(r, ranges.exons)
  dens[i] <- sum(q@width) / 50000
}
dens <- data.frame(x = s, dens = dens)

library(ggplot2)
library(cowplot)

p <- ggplot() + 
  #geom_rect(data = exons, aes(xmin = xmin, xmax = xmax), ymin = 0, ymax = 0.1, fill = "black") +
  geom_point(data = dat, aes(y=PI, x=BIN_START)) +
  geom_line(data = rec, aes(y=V2/300, x=V1), color = "red") +
  geom_line(data = tmrca, aes(y=AverageTmrca/6e6, x=Start), color = "blue") +
  geom_line(data = dens, aes(y=dens/50, x=x), color = "green")
p
ggsave(p, filename="Sim1.pdf", width = 8, height = 5)

p.rec <- ggplot(data = rec, aes(y=V2, x=V1), color = "red") + geom_line() + ylab("Recombination Rate") + xlab("Position")
p.tmr <- ggplot(data = tmrca, aes(y=AverageTmrca, x=Start), color = "blue") + geom_line() + xlab("Position")
p.dns <- ggplot(data = dens, aes(y=dens, x=x), color = "green") + geom_line() + ylab("Exonic site density") + xlab("Position")
p.taj <- ggplot(data = taj, aes(y=TajimaD, x=BIN_START), color = "orange") + geom_line() + ylab("Tajima's D") + xlab("Position")
p.snp <- ggplot(data = snp, aes(y=SNP_COUNT, x=BIN_START), color = "pink") + geom_line() + ylab("Nb. of variants per 50 kb") + xlab("Position")
p.m30 <- ggplot(data = m30, aes(y=SNP_COUNT, x=BIN_START), color = "pink") + geom_line() + ylab("Nb. of variants with MAF > 30% per 50 kb") + xlab("Position")
p.sin <- ggplot(data = sin, aes(y=SNP_COUNT, x=BIN_START), color = "purple") + geom_line() + ylab("Nb. of singletons per 50 kb") + xlab("Position")
p <- plot_grid(p.rec, p.tmr, p.dns, p.taj, p.m30, p.sin, nrow = 3, labels = "AUTO") 
ggsave(p, filename="Sim1.pdf", width = 12, height = 12)

#Debugging recapitation:

tmrca <- read.table("dm2L_bgs_rep1_w50000.tmrca_sample5.csv", header = T, sep = ",")
library(ggplot2)
library(reshape)
#dat <- melt(tmrca, measure.vars = c("AverageTmrcaOrig", "AverageTmrcaSimp", "AverageTmrcaRecap", "AverageTmrcaSimp5"), value.name = "TMRCA")
dat <- melt(tmrca, measure.vars = c("AverageTmrcaOrig", "AverageTmrcaRecap", "AverageTmrcaSimp5"), value.name = "TMRCA")

p <- ggplot(data = dat, aes(y=value, x=Start), color = "blue") + geom_line() + xlab("Position") + facet_wrap(~variable, ncol = 1)
p

ggsave(p, filename = "TMRCA_rep1nosimp.pdf", width = 8, height = 10)

#Old stuff:
dat <- read.table("Old/BGS_5.windowed.pi", header = TRUE)
#dat <- read.table("dm2L_bgs_old_w50000.windowed.pi", header = TRUE)
rec <- read.table("Comeron_tables/Comeron_100kb_chr2L.txt", header = F)

library(ggplot2)
p <- ggplot(dat, aes(y=PI, x=BIN_START)) + geom_point() +
  geom_line(dat = rec, aes(y=V2/10000, x=V1), color = "red")
p

```

