
library(plyr)
library(dplyr)

n <- 10 # number of replicates

avg_pairwise <- function(x) {
  
  pi <- 0
  
  if(any(x) == 1) {
    for(i in 1:n) {
      if(i < n) {
        for(j in (i+1):n) {
          if(x[i] != x[j]) {
            pi <- pi + 1
          }
        }
      }
    }
  }
  return(pi / choose(n, 2))
}

pb <- txtProgressBar(min = 0, max = n, style = 3)
for(r in 1:n) {
  
  setTxtProgressBar(pb, r)
  
  dat <- read.csv(gzfile(paste("rep_", r, "/rep_", r, ".gz", sep = "")), sep = " ")

  pi1 <- apply(dat, 1, avg_pairwise)
  
  pi_dat <- cbind.data.frame(pi1, 1:length(pi1))
  names(pi_dat) <- c("pi", "pos")
  pi_dat$bin.50kb <- (pi_dat$pos %/% (5e+4 + 1)) + 1
  pi_dat$bin.200kb <- (pi_dat$pos %/% (2e+5 + 1)) + 1
  pi_dat$bin.1Mb <- (pi_dat$pos %/% (1e+6 + 1)) + 1

  pi50k <- ddply(.data = pi_dat[,c(1, 3)], .variables = "bin.50kb", .fun = colMeans)
  pi50k$start <- seq(from = 0, length.out = nrow(pi50k), by = 5e+4)
  pi50k$end <- pi50k$start + 5e+4
  pi50k <- pi50k[,c(2:4, 1)]
  write.table(pi50k, paste("rep_", r, ".pi.50kb.bedgraph", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)

  pi200k <- ddply(.data = pi_dat[,c(1, 4)], .variables = "bin.200kb", .fun = colMeans)
  pi200k$start <- seq(from = 0, length.out = nrow(pi200k), by = 2e+5)
  pi200k$end <- pi200k$start + 2e+5
  pi200k <- pi200k[,c(2:4, 1)]
  write.table(pi200k, paste("rep_", r, ".pi.200kb.bedgraph", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)

  pi1M <- ddply(.data = pi_dat[,c(1, 5)], .variables = "bin.1Mb", .fun = colMeans)
  pi1M$start <- seq(from = 0, length.out = nrow(pi1M), by = 1e+6)
  pi1M$end <- pi1M$start + 1e+6
  pi1M <- pi1M[,c(2:4, 1)]
  write.table(pi1M, paste("rep_", r, ".pi.1Mb.bedgraph", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
}
close(pb)

