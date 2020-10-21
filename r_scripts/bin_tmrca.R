

library(plyr)
library(tidyverse)


args = commandArgs(trailingOnly = TRUE)

rep_idx <- as.character(args[1])



# gets nucleotide spans of each tree
arg <- readLines(paste("/rep_", rep_idx, "/rep_", rep_idx, ".newick", sep = ""))
tree.spans <- list()
for(i in 1:length(arg)) {
  tree.spans[i] <- data.matrix(strsplit(arg[[i]], "]", fixed = F))
}
nuc.spans <- character()
for(i in 1:length(tree.spans)) {
  nuc.spans[i] <- tree.spans[[i]][1]
}
tree.spans <- character()
for(i in 1:length(nuc.spans)) {
  tree.spans[i] <- data.matrix(strsplit(nuc.spans[[i]], "[", fixed = T))
}
for(i in 1:length(tree.spans)) {
  nuc.spans[i] <- tree.spans[[i]][2]
}
nuc.spans <- as.numeric(nuc.spans)
summary(nuc.spans)

sim.ARG <- read.tree(paste("/rep_", rep_idx, "/rep_", rep_idx, ".newick", sep = ""))

lst <- llply(sim.ARG, cophenetic.phylo, .progress = "text")

compute.mean.tmrca <- function(dist.mat) {
  row.sums <- apply(dist.mat, 1, sum)
  total.sum <- sum(row.sums)
  n <- ncol(dist.mat)
  return(total.sum / (n * n - n))
}

trees.tmrca <- llply(lst, compute.mean.tmrca, .progress = "text")
trees.tmrca <- as.numeric(trees.tmrca)

tree.seq <- as.data.frame(cbind(as.numeric(trees.tmrca), as.numeric(nuc.spans)))
row.names(tree.seq) <- 1:nrow(tree.seq)

# converts to single-nucleotide TMRCA landscape
tmrca.single.nuc <- list()

pb <- txtProgressBar(min = 0, max = nrow(tree.seq), style = 3)
for(i in 1:nrow(tree.seq)){
  setTxtProgressBar(pb, i)
  focal.height <- tree.seq[i, 1]
  focal.span <- tree.seq[i, 2]
  tmrca.single.nuc[[i]] <- rep(focal.height, focal.span)
}
close(pb)

tmrca.single.nuc <- do.call(c, tmrca.single.nuc)
tmrca.single.nuc <- as.data.frame(cbind(1:length(tmrca.single.nuc), tmrca.single.nuc))
names(tmrca.single.nuc) <- c("pos", "tmrca")

# 50kb
tmrca.single.nuc$bin <- ceiling(1:nrow(tmrca.single.nuc) / 50e+3)
sim.tmrca.50k <- ddply(.data = tmrca.single.nuc[-which(names(tmrca.single.nuc) == "pos")],
                       .variables = "bin", .fun = colMeans, .progress = "text")

write.table(sim.tmrca.50k, paste("rep_", rep_idx, "/sim.tmrca.50k.map", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)

# 200kb
tmrca.single.nuc$bin <- ceiling(1:nrow(tmrca.single.nuc) / 200e+3)
sim.tmrca.200k <- ddply(.data = tmrca.single.nuc[-which(names(tmrca.single.nuc) == "pos")],
                        .variables = "bin", .fun = colMeans, .progress = "text")

write.table(sim.tmrca.200k, paste("rep_", rep_idx, "/sim.tmrca.200k.map", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)

# 1Mb
tmrca.single.nuc$bin <- ceiling(1:nrow(tmrca.single.nuc) / 1e+6)
sim.tmrca.1M <- ddply(.data = tmrca.single.nuc[-which(names(tmrca.single.nuc) == "pos")],
                      .variables = "bin", .fun = colMeans, .progress = "text")

write.table(sim.tmrca.1M, paste("rep_", rep_idx, "/sim.tmrca.1M.map", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
