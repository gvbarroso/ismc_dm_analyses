
Generate a random mutation map:
```r
f <- "MutationMap.csv"
cat("Pos,RelRate\n", file = f)
pos <- 0
L <- 30e6
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
  window.ranges <- IRanges(start = seq(1, 30e6, by = window.size), width = window.size)
  ave <- numeric(length(window.ranges))
  for (i in 1:length(window.ranges)) {
    x <- restrict(rate.ranges, start(window.ranges[i]), end(window.ranges[i]))
    ave[i] <- weighted.mean(x@elementMetadata$Rate, w = width(x))
  }
  res <- data.frame(Start = start(window.ranges), Stop = end(window.ranges), Rate = ave)
  return(res)
}

mut <- read.table("MutationMap.csv", sep = ",", header = TRUE)
mut$Start <- mut$Pos + 1
mut$Stop <- c(mut$Pos[-1], 30e6)
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
rec <- data.frame(rates = c(0.1, 1.5, 0.1, 1.5, 0.1, 1.5) * 1e-8, start = c(0, 5, 10, 15, 20, 25) * 1e6)
rec.ranges <- IRanges(start = rec$start+1, width = 5e6, Rate = rec$rates)
rec.w50kb <- moving.average(rec.ranges, 50e3)
write.table(rec.w50kb, file = "RecombinationMap_50kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
rec.w200kb <- moving.average(rec.ranges, 200e3)
write.table(rec.w200kb, file = "RecombinationMap_200kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
rec.w1mb <- moving.average(rec.ranges, 1e6)
write.table(rec.w1mb, file = "RecombinationMap_1000kb.csv", sep = ",", row.names = FALSE, quote = FALSE)
```

