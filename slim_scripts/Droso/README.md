Get gene annotations:
```bash
mkdir GeneAnnotations
cd GeneAnnotations
wget http://ftp.ensembl.org/pub/release-103/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.103.gff3.gz
cd ..
```

We retrieve gene coordinates:
```r
library(ape)
gff <- read.gff("GeneAnnotations/Drosophila_melanogaster.BDGP6.32.103.gff3.gz")
table(gff$type)
exons <- subset(gff, type == "exon")
# Get chr 2L only:
exons.2L <- subset(exons, seqid == "2L")
# Collapse overlapping exons:
library(IRanges)
ranges <- IRanges(start = exons.2L$start, end = exons.2L$end)
ranges <- reduce(ranges)
exons.2L <- as.data.frame(ranges)
write.table(exons.2L[,c("start", "end")], file = "GeneAnnotations/exons_chr2L.csv", col.names = FALSE, row.names = FALSE, sep = ",")
```

Generate a random mutation map:
```r
f <- "MutationMap.csv"
cat("Pos,RelRate\n", file = f)
pos <- 0
L <- 23.51e6
while (pos < L) {
  relRate <- rgamma(1, shape = 2.5, rate = 2.5)
  cat(pos, ",", relRate, "\n", sep = "", file = f, append = TRUE)
  pos <- pos + rgeom(1, prob = 1e-5)
}
```

Get average mutation rates in windows:
```r
library(IRanges)
moving.average <- function(rate.ranges, window.size) {
  window.ranges <- IRanges(start = seq(1, 23.51e6, by = window.size), width = window.size)
  ave <- numeric(length(window.ranges))
  for (i in 1:length(window.ranges)) {
    x <- restrict(rate.ranges, start(window.ranges[i]), end(window.ranges[i]))
    ave[i] <- weighted.mean(x@elementMetadata$Rate, weights = width(x))
  }
  res <- data.frame(Start = start(window.ranges), Stop = end(window.ranges), Rate = ave)
  return(res)
}

mut <- read.table("MutationMap.csv", sep = ",", header = TRUE)
mut$Start <- mut$Pos + 1
mut$Stop <- c(mut$Pos[-1], 23.51e6)
mut.ranges <- IRanges(start = mut$Start, end = mut$Stop, Rate = mut$RelRate)

mut.w50kb <- moving.average(mut.ranges, 50e3)
write.table(mut.w50kb, file = "MutationMap_50kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
mut.w200kb <- moving.average(mut.ranges, 200e3)
write.table(mut.w200kb, file = "MutationMap_200kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
mut.w1mb <- moving.average(mut.ranges, 1e6)
write.table(mut.w1mb, file = "MutationMap_1000kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
```


Get average recombination rates in windows:
```r
rec <- read.table("Comeron_tables/Comeron_100kb_chr2L.txt", header = F)
rec.ranges <- IRanges(start = rec$V1, width = 100e3, Rate = rec$V2)

rec.w50kb <- moving.average(rec.ranges, 50e3)
write.table(rec.w50kb, file = "RecombinationMap_50kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
rec.w200kb <- moving.average(rec.ranges, 200e3)
write.table(rec.w200kb, file = "RecombinationMap_200kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
rec.w1mb <- moving.average(rec.ranges, 1e6)
write.table(rec.w1mb, file = "RecombinationMap_1000kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
```

