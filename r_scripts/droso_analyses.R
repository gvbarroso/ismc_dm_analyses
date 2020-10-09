# Created: 03/07/2019
# Last modified: 11/08/2020
# Author: Gustavo Barroso

library(ppcor)
library(MASS)
library(reshape2)
library(ggplot2)
library(cowplot)
library(lmtest)
library(nlme)
library(car)
library(plyr)
library(GenomicRanges)

setwd("~")
setwd("Data/iSMC/theta_paper/real_data/")

# R_2 table for plotting at the end
r2.chr.tab <- as.data.frame(matrix(ncol = 5, nrow = 3))
colnames(r2.chr.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Bin_size(kb)")

theme.blank <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

##################c####################
#
# Drosophila-like neutral simulations of 2L vs Real Data 2l
#
########################################

# building linear models

nreps <- 10
cnames <- character(length = nreps)
for(i in 1:nreps) { cnames[i] <- paste("rep_", i, sep = "") }

# 50 kb
r2.sim.tab.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.tab.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.tab.50kb) <- cnames

for(i in 1:nreps) {
  p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_", i, "/rs.pair_", i, ".", sep = "")
  pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
  pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
  tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
  tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
  rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
  theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)
  
  inf.lands.50k <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
  names(inf.lands.50k) <- c("diversity", "theta", "rho", "tmrca")
  
  m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands.50k)
  
  # type 2
  anova.diversity <- Anova(m.diversity)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.tab.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.sim.tab.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.tab.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.tab.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.tab.50kb$average <- rowMeans(r2.sim.tab.50kb)
r2.sim.tab.50kb <- transform(r2.sim.tab.50kb, sd=apply(r2.sim.tab.50kb, 1, sd, na.rm = TRUE))


# 200 kb
r2.sim.tab.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.tab.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.tab.200kb) <- cnames

for(i in 1:nreps)
{
  p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_", i, "/rs.pair_", i, ".", sep = "")
  pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
  pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
  tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
  tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
  rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
  theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)
  
  inf.lands.200k <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
  names(inf.lands.200k) <- c("diversity", "theta", "rho", "tmrca")
  
  m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands.200k)
  
  # type 2
  anova.diversity <- Anova(m.diversity)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.tab.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.sim.tab.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.tab.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.tab.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.tab.200kb$average <- rowMeans(r2.sim.tab.200kb)
r2.sim.tab.200kb <- transform(r2.sim.tab.200kb, sd=apply(r2.sim.tab.200kb, 1, sd, na.rm = TRUE))


# 1 Mb
r2.sim.tab.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.tab.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.tab.1Mb) <- cnames

for(i in 1:nreps)
{
  p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_", i, "/rs.pair_", i, ".", sep = "")
  pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
  pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
  tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
  tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
  rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
  theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)
  
  inf.lands.1M <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
  names(inf.lands.1M) <- c("diversity", "theta", "rho", "tmrca")
  
  m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands.1M)
  
  # type 2
  anova.diversity <- Anova(m.diversity)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.tab.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.sim.tab.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.tab.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.tab.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.tab.1Mb$average <- rowMeans(r2.sim.tab.1Mb)
r2.sim.tab.1Mb <- transform(r2.sim.tab.1Mb, sd=apply(r2.sim.tab.1Mb, 1, sd, na.rm = TRUE))

# merges
r2.sim.tab.avg <- rbind(r2.sim.tab.50kb$average, r2.sim.tab.200kb$average, r2.sim.tab.1Mb$average)
r2.sim.tab.avg <- as.data.frame(r2.sim.tab.avg)
names(r2.sim.tab.avg) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.tab.avg$bin.size <- c(50, 200, 1000)

r2.sim.tab.sd <- rbind(r2.sim.tab.50kb$sd, r2.sim.tab.200kb$sd, r2.sim.tab.1Mb$sd)
r2.sim.tab.sd <- as.data.frame(r2.sim.tab.sd)
names(r2.sim.tab.sd) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.tab.sd$bin.size <- c(50, 200, 1000)


# Real data 2L

r2.dm.tab <- as.data.frame(matrix(ncol = 5, nrow = 3))
colnames(r2.dm.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Bin_size(kb)")

# 50kb
# recombination landscapes
rho.dm.50kb <- read.table("dm_chr_maps/2L/dm_30x5x5.rho.50kb.bedgraph", header = T)

# diversity
diversity.dm.50kb <- read.table("dm_chr_maps/2L/dm_30x5x5.diversity.50kb.bedgraph", header = T)
diversity.dm.50kb$avg <- apply(diversity.dm.50kb[4:ncol(diversity.dm.50kb)], 1, mean)

# mutation landscapes
theta.dm.50kb <- read.table("dm_chr_maps/2L/dm_30x5x5.theta.50kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.50kb <- read.table("dm_chr_maps/2L/dm_30x5x5.TMRCA.50kb.bedgraph", header = T)
tmrca.dm.50kb$sample_mean <- apply(tmrca.dm.50kb[4:ncol(tmrca.dm.50kb)], 1, mean)

# missing data
missing.prop.50kb <- read.table("dm_chr_maps/2L/dm_30x5x5.missing.prop.50kb.bedgraph", header = T)
intersect.50kb <- apply(missing.prop.50kb[4:ncol(missing.prop.50kb)], 1, function(x) any(x > 0.1)) 

dm.lands.50kb <- as.data.frame(cbind(diversity.dm.50kb$avg,
                                     theta.dm.50kb$sample_mean,
                                     rho.dm.50kb$sample_mean,
                                     tmrca.dm.50kb$sample_mean))

names(dm.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
dm.lands.50kb$bin <- 1:nrow(dm.lands.50kb)
# filters based on missing data ( > 25% per window)
dm.lands.50kb <- dm.lands.50kb[which(intersect.50kb == F),]

# OLS
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.50kb)
bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)

# type 2 ANOVA
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.dm.tab[1,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3]) * 100,
                    anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, 50)

# 200kb
# recombination landscapes
rho.dm.200kb <- read.table("dm_chr_maps/2L/dm_30x5x5.rho.200kb.bedgraph", header = T)

# diversity
diversity.dm.200kb <- read.table("dm_chr_maps/2L/dm_30x5x5.diversity.200kb.bedgraph", header = T)
diversity.dm.200kb$avg <- apply(diversity.dm.200kb[4:ncol(diversity.dm.200kb)], 1, mean)

# mutation landscapes
theta.dm.200kb <- read.table("dm_chr_maps/2L/dm_30x5x5.theta.200kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.200kb <- read.table("dm_chr_maps/2L/dm_30x5x5.TMRCA.200kb.bedgraph", header = T)
tmrca.dm.200kb$sample_mean <- apply(tmrca.dm.200kb[4:ncol(tmrca.dm.200kb)], 1, mean)

# missing data
missing.prop.200kb <- read.table("dm_chr_maps/2L/dm_30x5x5.missing.prop.200kb.bedgraph", header = T)
intersect.200kb <- apply(missing.prop.200kb[4:ncol(missing.prop.200kb)], 1, function(x) any(x > 0.1)) 

dm.lands.200kb <- as.data.frame(cbind(diversity.dm.200kb$avg,
                                     theta.dm.200kb$sample_mean,
                                     rho.dm.200kb$sample_mean,
                                     tmrca.dm.200kb$sample_mean))

names(dm.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
dm.lands.200kb$bin <- 1:nrow(dm.lands.200kb)
# filters based on missing data ( > 25% per window)
dm.lands.200kb <- dm.lands.200kb[which(intersect.200kb == F),]

# OLS
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.200kb)
bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)

# type 2 ANOVA
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.dm.tab[2,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3]) * 100,
                    anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, 50)

# 1Mb
# recombination landscapes
rho.dm.1Mb <- read.table("dm_chr_maps/2L/dm_30x5x5.rho.1Mb.bedgraph", header = T)

# diversity
diversity.dm.1Mb <- read.table("dm_chr_maps/2L/dm_30x5x5.diversity.1Mb.bedgraph", header = T)
diversity.dm.1Mb$avg <- apply(diversity.dm.1Mb[4:ncol(diversity.dm.1Mb)], 1, mean)

# mutation landscapes
theta.dm.1Mb <- read.table("dm_chr_maps/2L/dm_30x5x5.theta.1Mb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.1Mb <- read.table("dm_chr_maps/2L/dm_30x5x5.TMRCA.1Mb.bedgraph", header = T)
tmrca.dm.1Mb$sample_mean <- apply(tmrca.dm.1Mb[4:ncol(tmrca.dm.1Mb)], 1, mean)

# missing data
missing.prop.1Mb <- read.table("dm_chr_maps/2L/dm_30x5x5.missing.prop.1Mb.bedgraph", header = T)
intersect.1Mb <- apply(missing.prop.1Mb[4:ncol(missing.prop.1Mb)], 1, function(x) any(x > 0.1)) 

dm.lands.1Mb <- as.data.frame(cbind(diversity.dm.1Mb$avg,
                                      theta.dm.1Mb$sample_mean,
                                      rho.dm.1Mb$sample_mean,
                                      tmrca.dm.1Mb$sample_mean))

names(dm.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
dm.lands.1Mb$bin <- 1:nrow(dm.lands.1Mb)
# filters based on missing data ( > 25% per window)
dm.lands.1Mb <- dm.lands.1Mb[which(intersect.1Mb == F),]

# OLS
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.1Mb)
bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)

# type 2 ANOVA
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)


r2.dm.tab[3,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3]) * 100,
                    anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, 1000)

# R2 Plot

r2.tab.2 <- as.data.frame(cbind(apply(r2.dm.tab, 2, as.numeric)))
names(r2.tab.2)[5] <- "bin.size"

r2.tab.comb <- rbind.data.frame(r2.tab.2, r2.sim.tab.avg)
r2.tab.comb$type <- c(rep("real", 3), rep("sim", 3))
r2.tab.comb$bin.size <- r2.tab.comb$bin.size

molten.r2 <- melt(r2.tab.comb, id.vars = c("bin.size", "type"))
r2.plot <- ggplot(data = molten.r2, aes(x = bin.size, y = value, colour = variable, fill = type))
r2.plot <- r2.plot + geom_line(data = molten.r2)
r2.plot <- r2.plot + geom_point(aes(shape = type, colour = variable), size = 7)
r2.plot <- r2.plot + scale_x_continuous(breaks = c(50, 200, 1000)) 
r2.plot <- r2.plot + scale_y_continuous(breaks = pretty_breaks())
r2.plot <- r2.plot + labs(title = NULL, x = "Bin Size (kb)", y = "Variance Explained (%)") + theme.blank
r2.plot <- r2.plot + theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))

leg <- get_legend(r2.plot + theme(legend.position="bottom"))

fig4ab <- plot_grid(cor.plot, r2.plot + no.legend, labels = "AUTO") # cor.plot is from script rs_samp_size.R

fig4 <- plot_grid(fig4ab, leg, nrow = 2, labels = NULL, rel_heights = c(1, 0.1), scale = c(1, 0.2))

cowplot::save_plot("Figure4.pdf", plot = fig4, base_width = 10, base_height = 6, device = "pdf")
 

######################################
#
# 30x5x5 --- 50kb windows
#
########################################

# Chr 2L

# recombination landscapes
rho.dm.50kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.rho.50kb.bedgraph", header = T)

# diversity
diversity.dm.50kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.diversity.50kb.bedgraph", header = T)
diversity.dm.50kb.2L$avg <- apply(diversity.dm.50kb.2L[4:ncol(diversity.dm.50kb.2L)], 1, mean)

# mutation landscapes
theta.dm.50kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.theta.50kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.50kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.TMRCA.50kb.bedgraph", header = T)
tmrca.dm.50kb.2L$sample_mean <- apply(tmrca.dm.50kb.2L[4:ncol(tmrca.dm.50kb.2L)], 1, mean)

# missing data
missing.prop.50kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.missing.prop.50kb.bedgraph", header = T)
intersect.50kb.2L <- apply(missing.prop.50kb.2L[4:ncol(missing.prop.50kb.2L)], 1, function(x) any(x > 0.1)) 

dm.lands.50kb.2L <- as.data.frame(cbind(diversity.dm.50kb.2L$chromStart,
                                        diversity.dm.50kb.2L$chromEnd,
                                        diversity.dm.50kb.2L$avg,
                                        theta.dm.50kb.2L$sample_mean,
                                        rho.dm.50kb.2L$sample_mean,
                                        tmrca.dm.50kb.2L$sample_mean))
names(dm.lands.50kb.2L) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.50kb.2L$bin <- 1:nrow(dm.lands.50kb.2L)

# filters based on missing data 
dm.lands.50kb.2L <- dm.lands.50kb.2L[which(intersect.50kb.2L == F),]

dm.lands.50kb.2L$chr <- 2

# Chr 2R

# recombination landscapes
rho.dm.50kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.rho.50kb.bedgraph", header = T)

# diversity
diversity.dm.50kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.diversity.50kb.bedgraph", header = T)
diversity.dm.50kb.2R$avg <- apply(diversity.dm.50kb.2R[4:ncol(diversity.dm.50kb.2R)], 1, mean)

# mutation landscapes
theta.dm.50kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.theta.50kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.50kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.TMRCA.50kb.bedgraph", header = T)
tmrca.dm.50kb.2R$sample_mean <- apply(tmrca.dm.50kb.2R[4:ncol(tmrca.dm.50kb.2R)], 1, mean)

# missing data
missing.prop.50kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.missing.prop.50kb.bedgraph", header = T)
intersect.50kb.2R <- apply(missing.prop.50kb.2R[4:ncol(missing.prop.50kb.2R)], 1, function(x) any(x > 0.1)) 

dm.lands.50kb.2R <- as.data.frame(cbind(diversity.dm.50kb.2R$chromStart,
                                        diversity.dm.50kb.2R$chromEnd,
                                        diversity.dm.50kb.2R$avg,
                                        theta.dm.50kb.2R$sample_mean,
                                        rho.dm.50kb.2R$sample_mean,
                                        tmrca.dm.50kb.2R$sample_mean))
names(dm.lands.50kb.2R) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.50kb.2R$bin <- 1:nrow(dm.lands.50kb.2R)

# filters based on missing data ( > 25% per window)
dm.lands.50kb.2R <- dm.lands.50kb.2R[which(intersect.50kb.2R == F),]

dm.lands.50kb.2R$chr <- 3

# Chr 3L

# recombination landscapes
rho.dm.50kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.rho.50kb.bedgraph", header = T)

# diversity
diversity.dm.50kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.diversity.50kb.bedgraph", header = T)
diversity.dm.50kb.3L$avg <- apply(diversity.dm.50kb.3L[4:ncol(diversity.dm.50kb.3L)], 1, mean)

# mutation landscapes
theta.dm.50kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.theta.50kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.50kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.TMRCA.50kb.bedgraph", header = T)
tmrca.dm.50kb.3L$sample_mean <- apply(tmrca.dm.50kb.3L[4:ncol(tmrca.dm.50kb.3L)], 1, mean)

# missing data
missing.prop.50kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.missing.prop.50kb.bedgraph", header = T)
intersect.50kb.3L <- apply(missing.prop.50kb.3L[4:ncol(missing.prop.50kb.3L)], 1, function(x) any(x > 0.1)) 

dm.lands.50kb.3L <- as.data.frame(cbind(diversity.dm.50kb.3L$chromStart,
                                        diversity.dm.50kb.3L$chromEnd,
                                        diversity.dm.50kb.3L$avg,
                                        theta.dm.50kb.3L$sample_mean,
                                        rho.dm.50kb.3L$sample_mean,
                                        tmrca.dm.50kb.3L$sample_mean))
names(dm.lands.50kb.3L) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.50kb.3L$bin <- 1:nrow(dm.lands.50kb.3L)

# filters based on missing data ( > 25% per window)
dm.lands.50kb.3L <- dm.lands.50kb.3L[which(intersect.50kb.3L == F),]

dm.lands.50kb.3L$chr <- 4

# Chr 3R

# recombination landscapes
rho.dm.50kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.rho.50kb.bedgraph", header = T)

# diversity
diversity.dm.50kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.diversity.50kb.bedgraph", header = T)
diversity.dm.50kb.3R$avg <- apply(diversity.dm.50kb.3R[4:ncol(diversity.dm.50kb.3R)], 1, mean)

# mutation landscapes
theta.dm.50kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.theta.50kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.50kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.TMRCA.50kb.bedgraph", header = T)
tmrca.dm.50kb.3R$sample_mean <- apply(tmrca.dm.50kb.3R[4:ncol(tmrca.dm.50kb.3R)], 1, mean)

# missing data
missing.prop.50kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.missing.prop.50kb.bedgraph", header = T)
intersect.50kb.3R <- apply(missing.prop.50kb[4:ncol(missing.prop.50kb.3R)], 1, function(x) any(x > 0.1)) 

dm.lands.50kb.3R <- as.data.frame(cbind(diversity.dm.50kb.3R$chromStart,
                                        diversity.dm.50kb.3R$chromEnd,
                                        diversity.dm.50kb.3R$avg,
                                        theta.dm.50kb.3R$sample_mean,
                                        rho.dm.50kb.3R$sample_mean,
                                        tmrca.dm.50kb.3R$sample_mean))
names(dm.lands.50kb.3R) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.50kb.3R$bin <- 1:nrow(dm.lands.50kb.3R)

# filters based on missing data ( > 25% per window)
dm.lands.50kb.3R <- dm.lands.50kb.3R[which(intersect.50kb.3R == F),]

dm.lands.50kb.3R$chr <- 5

# all together now
dm.lands.50kb <- rbind.data.frame(dm.lands.50kb.2L, dm.lands.50kb.2R, dm.lands.50kb.3L, dm.lands.50kb.3R)

write.table(dm.lands.50kb, "dm_chr_maps/dm.maps.50kb.tsv", quote = F, sep = "\t", row.names = F, col.names = T)


# Plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(dm.lands.50kb[c(3,7,8)], id.vars = c("bin", "chr"))
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.pi.50kb.pdf", device = "pdf")

molten.rho <- melt(dm.lands.50kb[c(5,7,8)], id.vars = c("bin", "chr"))
rho.map <- ggplot(data = molten.rho, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.rho.50kb.pdf", device = "pdf")

molten.theta <- melt(dm.lands.50kb[c(4,7,8)], id.vars = c("bin", "chr"))
theta.map <- ggplot(data = molten.theta, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.theta.50kb.pdf", device = "pdf")

molten.tmrca <- melt(dm.lands.50kb[c(6,7,8)], id.vars = c("bin", "chr"))
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.tau.50kb.pdf", device = "pdf")

# genome-wide correlations
cor.test(~rho + diversity, data = dm.lands.50kb, method = "spearman")
# 0.20 p-value = 2e-13
cor.test(~rho + tmrca, data = dm.lands.50kb, method = "spearman")
# 0.48 p-value < 2.2e-16
cor.test(~theta + tmrca, data = dm.lands.50kb, method = "spearman")
# 0.46 p-value < 2.2e-16

# Linear models
m.diversity <- lm(diversity ~ (theta + rho + tmrca) * chr, data = dm.lands.50kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.55, 0.95, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) #***
hist(resid(m.diversity.bc), nclass = 30)
hmctest(m.diversity.bc, nsim = 3000) # ***
dwtest(m.diversity.bc) # ***
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # ***
     
summary(m.diversity.bc)
# Coefficients:
#  Estimate Std. Error   t value Pr(>|t|)    
# (Intercept) -1.3400478  0.0011190 -1197.539   <2e-16 ***
# theta        3.0876273  0.0370939    83.238   <2e-16 ***
# rho          0.0225349  0.0133014     1.694   0.0905   
# tmrca        0.0334108  0.0015053    22.195   <2e-16 ***
# chr         -0.0001019  0.0002636    -0.386   0.6993    
# theta:chr    0.0165356  0.0096458     1.714   0.0867   
# rho:chr     -0.0082608  0.0036584    -2.258   0.0241 *  
# tmrca:chr    0.0001394  0.0003714     0.375   0.7075    

# Residual standard error: 0.0009126 on 1288 degrees of freedom
# Multiple R-squared:  0.9929,	Adjusted R-squared:  0.9929 
# F-statistic: 2.584e+04 on 7 and 1288 DF,  p-value: < 2.2e-16

# type 2 ANOVA
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.tab[1,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3]) * 100,
                 anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, 50)

# model w/o chr to estimate VIF
m.no.interaction <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.50kb)
summary(m.no.interaction)
vif(m.no.interaction)
# theta      rho    tmrca 
# 1.29 1.57 1.93

# Because of auto-correlation we compute p-values for the variables using a GLS
dm.lands.50kb$idx <- 1:nrow(dm.lands.50kb)
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = dm.lands.50kb, corr = corAR1(0, ~idx))

summary(g.diversity)
# Parameter estimate(s):
# Phi 
# 0.1717399 

# Coefficients:
#  Value   Std.Error   t-value p-value
# (Intercept) -0.0086368 0.000105886 -81.56718  0.0000
# theta        0.9812938 0.004292210 228.62201  0.0000
# rho          0.0015014 0.001546371   0.97092  0.3318
# tmrca        0.0092748 0.000148631  62.40164  0.0000


# Linear model without TMRCA --> rho becomes significant
m.diversity.no.tau <- lm(diversity ~ theta + rho, data = dm.lands.50kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity.no.tau, lambda = seq(0.9, 1.3, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.no.tau.bc <- update(m.diversity.no.tau, (diversity^l -1)/l~.)
plot(m.diversity.no.tau.bc, which = 2)

shapiro.test(resid(m.diversity.no.tau.bc)) #***
hist(resid(m.diversity.no.tau.bc), nclass = 30)
hmctest(m.diversity.no.tau.bc, nsim = 3000) # ***
dwtest(m.diversity.no.tau.bc) # ***
Box.test(resid(m.diversity.no.tau.bc)[order(predict(m.diversity.no.tau.bc))], type = "Ljung-Box") # NS

summary(m.diversity.no.tau.bc)
#Coefficients:
#  Estimate Std. Error   t value Pr(>|t|)    
# (Intercept) -9.063e-01  5.912e-05 -15328.16   <2e-16 ***
# theta        6.508e-01  4.048e-03    160.78   <2e-16 ***
# rho          3.409e-02  1.433e-03     23.79   <2e-16 ***

# Residual standard error: 0.0004113 on 1293 degrees of freedom
# Multiple R-squared:  0.9556,	Adjusted R-squared:  0.9556 
# F-statistic: 1.393e+04 on 2 and 1293 DF,  p-value: < 2.2e-16

######################################
#
# 30x5x5 --- 200kb windows
#
########################################

# Chr 2L

# recombination landscapes
rho.dm.200kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.rho.200kb.bedgraph", header = T)

# diversity
diversity.dm.200kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.diversity.200kb.bedgraph", header = T)
diversity.dm.200kb.2L$avg <- apply(diversity.dm.200kb.2L[4:ncol(diversity.dm.200kb.2L)], 1, mean)

# mutation landscapes
theta.dm.200kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.theta.200kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.200kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.TMRCA.200kb.bedgraph", header = T)
tmrca.dm.200kb.2L$sample_mean <- apply(tmrca.dm.200kb.2L[4:ncol(tmrca.dm.200kb.2L)], 1, mean)

# missing data
missing.prop.200kb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.missing.prop.200kb.bedgraph", header = T)
intersect.200kb.2L <- apply(missing.prop.200kb.2L[4:ncol(missing.prop.200kb.2L)], 1, function(x) any(x > 0.1)) 

dm.lands.200kb.2L <- as.data.frame(cbind(diversity.dm.200kb.2L$chromStart,
                                        diversity.dm.200kb.2L$chromEnd,
                                        diversity.dm.200kb.2L$avg,
                                        theta.dm.200kb.2L$sample_mean,
                                        rho.dm.200kb.2L$sample_mean,
                                        tmrca.dm.200kb.2L$sample_mean))
names(dm.lands.200kb.2L) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.200kb.2L$bin <- 1:nrow(dm.lands.200kb.2L)

# filters based on missing data 
dm.lands.200kb.2L <- dm.lands.200kb.2L[which(intersect.200kb.2L == F),]

dm.lands.200kb.2L$chr <- 2

# Chr 2R

# recombination landscapes
rho.dm.200kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.rho.200kb.bedgraph", header = T)

# diversity
diversity.dm.200kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.diversity.200kb.bedgraph", header = T)
diversity.dm.200kb.2R$avg <- apply(diversity.dm.200kb.2R[4:ncol(diversity.dm.200kb.2R)], 1, mean)

# mutation landscapes
theta.dm.200kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.theta.200kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.200kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.TMRCA.200kb.bedgraph", header = T)
tmrca.dm.200kb.2R$sample_mean <- apply(tmrca.dm.200kb.2R[4:ncol(tmrca.dm.200kb.2R)], 1, mean)

# missing data
missing.prop.200kb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.missing.prop.200kb.bedgraph", header = T)
intersect.200kb.2R <- apply(missing.prop.200kb.2R[4:ncol(missing.prop.200kb.2R)], 1, function(x) any(x > 0.1)) 

dm.lands.200kb.2R <- as.data.frame(cbind(diversity.dm.200kb.2R$chromStart,
                                        diversity.dm.200kb.2R$chromEnd,
                                        diversity.dm.200kb.2R$avg,
                                        theta.dm.200kb.2R$sample_mean,
                                        rho.dm.200kb.2R$sample_mean,
                                        tmrca.dm.200kb.2R$sample_mean))
names(dm.lands.200kb.2R) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.200kb.2R$bin <- 1:nrow(dm.lands.200kb.2R)

# filters based on missing data ( > 25% per window)
dm.lands.200kb.2R <- dm.lands.200kb.2R[which(intersect.200kb.2R == F),]

dm.lands.200kb.2R$chr <- 3

# Chr 3L

# recombination landscapes
rho.dm.200kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.rho.200kb.bedgraph", header = T)

# diversity
diversity.dm.200kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.diversity.200kb.bedgraph", header = T)
diversity.dm.200kb.3L$avg <- apply(diversity.dm.200kb.3L[4:ncol(diversity.dm.200kb.3L)], 1, mean)

# mutation landscapes
theta.dm.200kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.theta.200kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.200kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.TMRCA.200kb.bedgraph", header = T)
tmrca.dm.200kb.3L$sample_mean <- apply(tmrca.dm.200kb.3L[4:ncol(tmrca.dm.200kb.3L)], 1, mean)

# missing data
missing.prop.200kb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.missing.prop.200kb.bedgraph", header = T)
intersect.200kb.3L <- apply(missing.prop.200kb.3L[4:ncol(missing.prop.200kb.3L)], 1, function(x) any(x > 0.1)) 

dm.lands.200kb.3L <- as.data.frame(cbind(diversity.dm.200kb.3L$chromStart,
                                        diversity.dm.200kb.3L$chromEnd,
                                        diversity.dm.200kb.3L$avg,
                                        theta.dm.200kb.3L$sample_mean,
                                        rho.dm.200kb.3L$sample_mean,
                                        tmrca.dm.200kb.3L$sample_mean))
names(dm.lands.200kb.3L) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.200kb.3L$bin <- 1:nrow(dm.lands.200kb.3L)

# filters based on missing data ( > 25% per window)
dm.lands.200kb.3L <- dm.lands.200kb.3L[which(intersect.200kb.3L == F),]

dm.lands.200kb.3L$chr <- 4

# Chr 3R

# recombination landscapes
rho.dm.200kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.rho.200kb.bedgraph", header = T)

# diversity
diversity.dm.200kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.diversity.200kb.bedgraph", header = T)
diversity.dm.200kb.3R$avg <- apply(diversity.dm.200kb.3R[4:ncol(diversity.dm.200kb.3R)], 1, mean)

# mutation landscapes
theta.dm.200kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.theta.200kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.200kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.TMRCA.200kb.bedgraph", header = T)
tmrca.dm.200kb.3R$sample_mean <- apply(tmrca.dm.200kb.3R[4:ncol(tmrca.dm.200kb.3R)], 1, mean)

# missing data
missing.prop.200kb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.missing.prop.200kb.bedgraph", header = T)
intersect.200kb.3R <- apply(missing.prop.200kb[4:ncol(missing.prop.200kb.3R)], 1, function(x) any(x > 0.1)) 

dm.lands.200kb.3R <- as.data.frame(cbind(diversity.dm.200kb.3R$chromStart,
                                        diversity.dm.200kb.3R$chromEnd,
                                        diversity.dm.200kb.3R$avg,
                                        theta.dm.200kb.3R$sample_mean,
                                        rho.dm.200kb.3R$sample_mean,
                                        tmrca.dm.200kb.3R$sample_mean))
names(dm.lands.200kb.3R) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.200kb.3R$bin <- 1:nrow(dm.lands.200kb.3R)

# filters based on missing data ( > 25% per window)
dm.lands.200kb.3R <- dm.lands.200kb.3R[which(intersect.200kb.3R == F),]

dm.lands.200kb.3R$chr <- 5

# all together now
dm.lands.200kb <- rbind.data.frame(dm.lands.200kb.2L, dm.lands.200kb.2R, dm.lands.200kb.3L, dm.lands.200kb.3R)

write.table(dm.lands.200kb, "dm_chr_maps/dm.maps.200kb.tsv", quote = F, sep = "\t", row.names = F, col.names = T)


# Plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(dm.lands.200kb[c(3,7,8)], id.vars = c("bin", "chr"))
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.pi.200kb.pdf", device = "pdf")

molten.rho <- melt(dm.lands.200kb[c(5,7,8)], id.vars = c("bin", "chr"))
rho.map <- ggplot(data = molten.rho, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.rho.200kb.pdf", device = "pdf")

molten.theta <- melt(dm.lands.200kb[c(4,7,8)], id.vars = c("bin", "chr"))
theta.map <- ggplot(data = molten.theta, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.theta.200kb.pdf", device = "pdf")

molten.tmrca <- melt(dm.lands.200kb[c(6,7,8)], id.vars = c("bin", "chr"))
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.tau.200kb.pdf", device = "pdf")

# genome-wide correlations
cor.test(~rho + diversity, data = dm.lands.200kb, method = "spearman")
# 0.15 p-value = 0.0025
cor.test(~rho + tmrca, data = dm.lands.200kb, method = "spearman")
# 0.45 p-value < 2.2e-16
cor.test(~theta + tmrca, data = dm.lands.200kb, method = "spearman")
# 0.51 p-value < 2.2e-16

# Linear models
m.diversity <- lm(diversity ~ (theta + rho + tmrca) * chr, data = dm.lands.200kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.6, 1.1, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # NS
hist(resid(m.diversity.bc), nclass = 30)
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # ***
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # ***

summary(m.diversity.bc)
# Coefficients:
# Estimate Std. Error   t value Pr(>|t|)    
# (Intercept) -1.198e+00  1.065e-03 -1124.544   <2e-16 ***
# theta        2.097e+00  2.726e-02    76.919   <2e-16 ***
# rho          2.205e-02  1.193e-02     1.848   0.0655   
# tmrca        1.845e-02  1.449e-03    12.740   <2e-16 ***
# chr         -4.206e-05  2.379e-04    -0.177   0.8598    
# theta:chr    1.530e-03  6.817e-03     0.224   0.8225    
# rho:chr     -4.551e-03  3.051e-03    -1.492   0.1366    
# tmrca:chr    1.560e-04  3.354e-04     0.465   0.6420    

# Residual standard error: 0.0003042 on 354 degrees of freedom
# Multiple R-squared:  0.9978,	Adjusted R-squared:  0.9977 
# F-statistic: 2.277e+04 on 7 and 354 DF,  p-value: < 2.2e-16

# type 2 ANOVA
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.tab[2,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3]) * 100,
                anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, 200)

# model w/o chr to estimate VIF
m.no.interaction <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.200kb)
summary(m.no.interaction)
vif(m.no.interaction)
# theta      rho    tmrca 
# 1.51 1.61 2.24

# Because of auto-correlation we compute p-values for the variables using a GLS
dm.lands.200kb$idx <- 1:nrow(dm.lands.200kb)
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = dm.lands.200kb, corr = corAR1(0, ~idx))

summary(g.diversity)
# Parameter estimate(s):
# Phi 
# 0.297378 

# Coefficients:
# Value   Std.Error   t-value p-value
# (Intercept) -0.0076045 0.000147761 -51.46460  0.0000
# theta        0.9929796 0.005014208 198.03320  0.0000
# rho          0.0053508 0.002126439   2.51632  0.0123
# tmrca        0.0079198 0.000210546  37.61566  0.0000


# Linear model without TMRCA --> rho becomes significant
m.diversity.no.tau <- lm(diversity ~ theta + rho, data = dm.lands.200kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity.no.tau, lambda = seq(1, 1.4, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.no.tau.bc <- update(m.diversity.no.tau, (diversity^l -1)/l~.)
plot(m.diversity.no.tau.bc, which = 2)

shapiro.test(resid(m.diversity.no.tau.bc)) #***
hist(resid(m.diversity.no.tau.bc), nclass = 30)
hmctest(m.diversity.no.tau.bc, nsim = 3000) # **
dwtest(m.diversity.no.tau.bc) # ***
Box.test(resid(m.diversity.no.tau.bc)[order(predict(m.diversity.no.tau.bc))], type = "Ljung-Box") # NS

summary(m.diversity.no.tau.bc)
# Coefficients:
# Estimate Std. Error   t value Pr(>|t|)    
# (Intercept) -8.395e-01  5.147e-05 -16311.40   <2e-16 ***
# theta        4.222e-01  2.953e-03    142.97   <2e-16 ***
# rho          1.897e-02  1.329e-03     14.27   <2e-16 ***
  
# Residual standard error: 0.0001468 on 359 degrees of freedom
# Multiple R-squared:  0.9836,	Adjusted R-squared:  0.9835 
# F-statistic: 1.079e+04 on 2 and 359 DF,  p-value: < 2.2e-16

######################################
#
# 30x5x5 --- 1Mb windows
#
########################################

# Chr 2L

# recombination landscapes
rho.dm.1Mb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.rho.1Mb.bedgraph", header = T)

# diversity
diversity.dm.1Mb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.diversity.1Mb.bedgraph", header = T)
diversity.dm.1Mb.2L$avg <- apply(diversity.dm.1Mb.2L[4:ncol(diversity.dm.1Mb.2L)], 1, mean)

# mutation landscapes
theta.dm.1Mb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.theta.1Mb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.1Mb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.TMRCA.1Mb.bedgraph", header = T)
tmrca.dm.1Mb.2L$sample_mean <- apply(tmrca.dm.1Mb.2L[4:ncol(tmrca.dm.1Mb.2L)], 1, mean)

# missing data
missing.prop.1Mb.2L <- read.table("dm_chr_maps/2L/dm_30x5x5.missing.prop.1Mb.bedgraph", header = T)
intersect.1Mb.2L <- apply(missing.prop.1Mb.2L[4:ncol(missing.prop.1Mb.2L)], 1, function(x) any(x > 0.1)) 

dm.lands.1Mb.2L <- as.data.frame(cbind(diversity.dm.1Mb.2L$chromStart,
                                       diversity.dm.1Mb.2L$chromEnd,
                                       diversity.dm.1Mb.2L$avg,
                                       theta.dm.1Mb.2L$sample_mean,
                                         rho.dm.1Mb.2L$sample_mean,
                                         tmrca.dm.1Mb.2L$sample_mean))
names(dm.lands.1Mb.2L) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.1Mb.2L$bin <- 1:nrow(dm.lands.1Mb.2L)

# filters based on missing data ( > 25% per window)
dm.lands.1Mb.2L <- dm.lands.1Mb.2L[which(intersect.1Mb.2L == F),]

dm.lands.1Mb.2L$chr <- 2

# Chr 2R

# recombination landscapes
rho.dm.1Mb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.rho.1Mb.bedgraph", header = T)

# diversity
diversity.dm.1Mb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.diversity.1Mb.bedgraph", header = T)
diversity.dm.1Mb.2R$avg <- apply(diversity.dm.1Mb.2R[4:ncol(diversity.dm.1Mb.2R)], 1, mean)

# mutation landscapes
theta.dm.1Mb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.theta.1Mb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.1Mb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.TMRCA.1Mb.bedgraph", header = T)
tmrca.dm.1Mb.2R$sample_mean <- apply(tmrca.dm.1Mb.2R[4:ncol(tmrca.dm.1Mb.2R)], 1, mean)

# missing data
missing.prop.1Mb.2R <- read.table("dm_chr_maps/2R/dm_30x5x5.missing.prop.1Mb.bedgraph", header = T)
intersect.1Mb.2R <- apply(missing.prop.1Mb.2R[4:ncol(missing.prop.1Mb.2R)], 1, function(x) any(x > 0.1)) 

dm.lands.1Mb.2R <- as.data.frame(cbind(diversity.dm.1Mb.2R$chromStart,
                                         diversity.dm.1Mb.2R$chromEnd,
                                         diversity.dm.1Mb.2R$avg,
                                         theta.dm.1Mb.2R$sample_mean,
                                         rho.dm.1Mb.2R$sample_mean,
                                         tmrca.dm.1Mb.2R$sample_mean))
names(dm.lands.1Mb.2R) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.1Mb.2R$bin <- 1:nrow(dm.lands.1Mb.2R)

# filters based on missing data ( > 25% per window)
dm.lands.1Mb.2R <- dm.lands.1Mb.2R[which(intersect.1Mb.2R == F),]

dm.lands.1Mb.2R$chr <- 3

# Chr 3L

# recombination landscapes
rho.dm.1Mb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.rho.1Mb.bedgraph", header = T)

# diversity
diversity.dm.1Mb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.diversity.1Mb.bedgraph", header = T)
diversity.dm.1Mb.3L$avg <- apply(diversity.dm.1Mb.3L[4:ncol(diversity.dm.1Mb.3L)], 1, mean)

# mutation landscapes
theta.dm.1Mb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.theta.1Mb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.1Mb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.TMRCA.1Mb.bedgraph", header = T)
tmrca.dm.1Mb.3L$sample_mean <- apply(tmrca.dm.1Mb.3L[4:ncol(tmrca.dm.1Mb.3L)], 1, mean)

# missing data
missing.prop.1Mb.3L <- read.table("dm_chr_maps/3L/dm_30x5x5.missing.prop.1Mb.bedgraph", header = T)
intersect.1Mb.3L <- apply(missing.prop.1Mb.3L[4:ncol(missing.prop.1Mb.3L)], 1, function(x) any(x > 0.1)) 

dm.lands.1Mb.3L <- as.data.frame(cbind(diversity.dm.1Mb.3L$chromStart,
                                         diversity.dm.1Mb.3L$chromEnd,
                                         diversity.dm.1Mb.3L$avg,
                                         theta.dm.1Mb.3L$sample_mean,
                                         rho.dm.1Mb.3L$sample_mean,
                                         tmrca.dm.1Mb.3L$sample_mean))
names(dm.lands.1Mb.3L) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.1Mb.3L$bin <- 1:nrow(dm.lands.1Mb.3L)

# filters based on missing data 
dm.lands.1Mb.3L <- dm.lands.1Mb.3L[which(intersect.1Mb.3L == F),]

dm.lands.1Mb.3L$chr <- 4

# Chr 3R

# recombination landscapes
rho.dm.1Mb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.rho.1Mb.bedgraph", header = T)

# diversity
diversity.dm.1Mb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.diversity.1Mb.bedgraph", header = T)
diversity.dm.1Mb.3R$avg <- apply(diversity.dm.1Mb.3R[4:ncol(diversity.dm.1Mb.3R)], 1, mean)

# mutation landscapes
theta.dm.1Mb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.theta.1Mb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.1Mb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.TMRCA.1Mb.bedgraph", header = T)
tmrca.dm.1Mb.3R$sample_mean <- apply(tmrca.dm.1Mb.3R[4:ncol(tmrca.dm.1Mb.3R)], 1, mean)

# missing data
missing.prop.1Mb.3R <- read.table("dm_chr_maps/3R/dm_30x5x5.missing.prop.1Mb.bedgraph", header = T)
intersect.1Mb.3R <- apply(missing.prop.1Mb[4:ncol(missing.prop.1Mb.3R)], 1, function(x) any(x > 0.1)) 

dm.lands.1Mb.3R <- as.data.frame(cbind(diversity.dm.1Mb.3R$chromStart,
                                         diversity.dm.1Mb.3R$chromEnd,
                                         diversity.dm.1Mb.3R$avg,
                                         theta.dm.1Mb.3R$sample_mean,
                                         rho.dm.1Mb.3R$sample_mean,
                                         tmrca.dm.1Mb.3R$sample_mean))
names(dm.lands.1Mb.3R) <- c("start", "end", "diversity", "theta", "rho", "tmrca")
dm.lands.1Mb.3R$bin <- 1:nrow(dm.lands.1Mb.3R)

# filters based on missing data ( > 25% per window)
dm.lands.1Mb.3R <- dm.lands.1Mb.3R[which(intersect.1Mb.3R == F),]

dm.lands.1Mb.3R$chr <- 5

# all together now
dm.lands.1Mb <- rbind.data.frame(dm.lands.1Mb.2L, dm.lands.1Mb.2R, dm.lands.1Mb.3L, dm.lands.1Mb.3R)

write.table(dm.lands.1Mb, "dm_chr_maps/dm.maps.1Mb.tsv", quote = F, sep = "\t", row.names = F, col.names = T)


# Plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(dm.lands.1Mb[c(3,7,8)], id.vars = c("bin", "chr"))
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.pi.1Mb.pdf", device = "pdf")

molten.rho <- melt(dm.lands.1Mb[c(5,7,8)], id.vars = c("bin", "chr"))
rho.map <- ggplot(data = molten.rho, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.rho.1Mb.pdf", device = "pdf")

molten.theta <- melt(dm.lands.1Mb[c(4,7,8)], id.vars = c("bin", "chr"))
theta.map <- ggplot(data = molten.theta, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.theta.1Mb.pdf", device = "pdf")

molten.tmrca <- melt(dm.lands.1Mb[c(6,7,8)], id.vars = c("bin", "chr"))
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 50, y = value)) + facet_wrap(~chr)
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
ggsave("dm_chr_maps/dm.tau.1Mb.pdf", device = "pdf")

# genome-wide correlations
cor.test(~rho + diversity, data = dm.lands.1Mb, method = "spearman")
# 0.20 p-value = 0.07
cor.test(~rho + tmrca, data = dm.lands.1Mb, method = "spearman")
# 0.48 p-value < 2.2e-16
cor.test(~theta + tmrca, data = dm.lands.1Mb, method = "spearman")
# 0.66 p-value < 2.2e-16

# Linear models
m.diversity <- lm(diversity ~ (theta + rho + tmrca) * chr, data = dm.lands.1Mb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.6, 1.1, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # NS
hist(resid(m.diversity.bc), nclass = 30)
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # NS
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# Coefficients:
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept) -1.1108861  0.0019736 -562.881  < 2e-16 ***
# theta        1.5266659  0.0420695   36.289  < 2e-16 ***
# rho          0.0240532  0.0152588    1.576    0.119    
# tmrca        0.0116707  0.0025143    4.642 1.37e-05 ***
# chr         -0.0002336  0.0004183   -0.558    0.578    
# theta:chr    0.0067786  0.0096641    0.701    0.485    
# rho:chr     -0.0036485  0.0039469   -0.924    0.358    
# tmrca:chr    0.0002839  0.0005504    0.516    0.607    

# Residual standard error: 0.000148 on 78 degrees of freedom
# Multiple R-squared:  0.9989,	Adjusted R-squared:  0.9988 
# F-statistic: 1.022e+04 on 7 and 78 DF,  p-value: < 2.2e-16

# type 2 ANOVA
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.tab[3,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3]) * 100,
                 anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, 1000)

# model w/o chr to estimate VIF
m.no.interaction <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.1Mb)
summary(m.no.interaction)
vif(m.no.interaction)
# theta      rho    tmrca 
# 1.88 1.61 2.69

# Because of auto-correlation we compute p-values for the variables using a GLS
dm.lands.1Mb$idx <- 1:nrow(dm.lands.1Mb)
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = dm.lands.1Mb, corr = corAR1(0, ~idx))


# Linear model without TMRCA --> rho becomes significant
m.diversity.no.tau <- lm(diversity ~ theta + rho, data = dm.lands.1Mb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity.no.tau, lambda = seq(1, 1.4, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.no.tau.bc <- update(m.diversity.no.tau, (diversity^l -1)/l~.)
plot(m.diversity.no.tau.bc, which = 2)

shapiro.test(resid(m.diversity.no.tau.bc)) # NS
hist(resid(m.diversity.no.tau.bc), nclass = 30)
hmctest(m.diversity.no.tau.bc, nsim = 3000) # NS
dwtest(m.diversity.no.tau.bc) # ***
Box.test(resid(m.diversity.no.tau.bc)[order(predict(m.diversity.no.tau.bc))], type = "Ljung-Box") # NS

summary(m.diversity.no.tau.bc)
# Coefficients:
#  Estimate Std. Error    t value Pr(>|t|)    
# (Intercept) -7.899e-01  6.069e-05 -13014.329  < 2e-16 ***
# theta        2.961e-01  3.265e-03     90.675  < 2e-16 ***
# rho          1.248e-02  1.649e-03      7.573 4.58e-11 ***

# Residual standard error: 7.02e-05 on 83 degrees of freedom
# Multiple R-squared:  0.9908,	Adjusted R-squared:  0.9906 
# F-statistic:  4482 on 2 and 83 DF,  p-value: < 2.2e-16

########################################
#
# Divergence --- 50kb maps for ++ resolution
#
########################################

# divergence data from D. melanogaster and D. yakuba
divergence.2L.5kb <- read.table("~/Data/iSMC/theta_paper/real_data/dm_misc/Droso2L_divergence.statistics5kb.csv", header = T)
divergence.2R.5kb <- read.table("~/Data/iSMC/theta_paper/real_data/dm_misc/Droso2R_divergence.statistics5kb.csv", header = T)
divergence.3L.5kb <- read.table("~/Data/iSMC/theta_paper/real_data/dm_misc/Droso3L_divergence.statistics5kb.csv", header = T)
divergence.3R.5kb <- read.table("~/Data/iSMC/theta_paper/real_data/dm_misc/Droso3R_divergence.statistics5kb.csv", header = T)

divergence <- rbind.data.frame(divergence.2L.5kb, divergence.2R.5kb, divergence.3L.5kb, divergence.3R.5kb)
divergence <- divergence[,c(1:3, 6)] 
divergence$Chr <- as.character(divergence$Chr)
  
# reorder maps to use GR
dm.maps.50kb <- dm.lands.50kb[,c(8,1:6)]
dm.maps.50kb$chr[which(dm.maps.50kb$chr == 2)] <- "2L"
dm.maps.50kb$chr[which(dm.maps.50kb$chr == 3)] <- "2R"
dm.maps.50kb$chr[which(dm.maps.50kb$chr == 4)] <- "3L"
dm.maps.50kb$chr[which(dm.maps.50kb$chr == 5)] <- "3R"

# converts objects to GenomicRanges
dm.gr <- makeGRangesFromDataFrame(dm.maps.50kb)
values(dm.gr) <- dm.maps.50kb[,(4:7)]
divergence.gr <- makeGRangesFromDataFrame(divergence)
values(divergence.gr) <- DataFrame(score = divergence$MLModelFit.BrLen0)

hits <- findOverlaps(query = divergence.gr, subject = dm.gr, type = "within")
ranges(divergence.gr)[queryHits(hits)] = ranges(dm.gr)[subjectHits(hits)]

lands.gr.df <- as.data.frame(dm.gr)
divergence.gr.df <- as.data.frame(divergence.gr)
# deletes non-matching windows
divergence.gr.df <- divergence.gr.df[-which(((divergence.gr.df$width - 1) %% 50000) != 0),]

# compute mean divergence within 50kb windows
dummy.tbl <- divergence.gr.df[, -c(4, 5)]
dummy.tbl$seqnames <- as.character(dummy.tbl$seqnames)
dummy.tbl$seqnames[which(dummy.tbl$seqnames == "2L")] <- 2
dummy.tbl$seqnames[which(dummy.tbl$seqnames == "2R")] <- 3
dummy.tbl$seqnames[which(dummy.tbl$seqnames == "3L")] <- 4
dummy.tbl$seqnames[which(dummy.tbl$seqnames == "3R")] <- 5
dummy.tbl$seqnames <- as.numeric(dummy.tbl$seqnames)

tmp <- ddply(.data = dummy.tbl, .variables = c("seqnames", "start"), .fun = colMeans, .progress = "text")

divergence.gr.df.2L <- divergence.gr.df[divergence.gr.df$seqnames == "2L",] 
divergence.gr.df.2L <- divergence.gr.df.2L[!duplicated(divergence.gr.df.2L$start),]
divergence.gr.df.2R <- divergence.gr.df[divergence.gr.df$seqnames == "2R",] 
divergence.gr.df.2R <- divergence.gr.df.2R[!duplicated(divergence.gr.df.2R$start),]
divergence.gr.df.3L <- divergence.gr.df[divergence.gr.df$seqnames == "3L",] 
divergence.gr.df.3L <- divergence.gr.df.3L[!duplicated(divergence.gr.df.3L$start),]
divergence.gr.df.3R <- divergence.gr.df[divergence.gr.df$seqnames == "3R",] 
divergence.gr.df.3R <- divergence.gr.df.3R[!duplicated(divergence.gr.df.3R$start),]

divergence.gr.df.chr <- rbind.data.frame(divergence.gr.df.2L, divergence.gr.df.2R, divergence.gr.df.3L, divergence.gr.df.3R)
divergence.gr.df.chr$score <- tmp$score

lands.gr.df.2L <- lands.gr.df[which(lands.gr.df$seqnames == "2L"),]
lands.gr.df.2L <- lands.gr.df.2L[which(lands.gr.df.2L$start %in% divergence.gr.df.2L$start),]
lands.gr.df.2R <- lands.gr.df[which(lands.gr.df$seqnames == "2R"),]
lands.gr.df.2R <- lands.gr.df.2R[which(lands.gr.df.2R$start %in% divergence.gr.df.2R$start),]
lands.gr.df.3L <- lands.gr.df[which(lands.gr.df$seqnames == "3L"),]
lands.gr.df.3L <- lands.gr.df.3L[which(lands.gr.df.3L$start %in% divergence.gr.df.3L$start),]
lands.gr.df.3R <- lands.gr.df[which(lands.gr.df$seqnames == "3R"),]
lands.gr.df.3R <- lands.gr.df.3R[which(lands.gr.df.3R$start %in% divergence.gr.df.3R$start),]

lands.gr.df.chr <- rbind.data.frame(lands.gr.df.2L, lands.gr.df.2R, lands.gr.df.3L, lands.gr.df.3R)

lands.divergence.dm <- cbind.data.frame(lands.gr.df.chr[,-(which(names(lands.gr.df.chr) == "strand"))], divergence.gr.df.chr$score)
names(lands.divergence.dm)[ncol(lands.divergence.dm)] <- "divergence"

write.table(lands.divergence.dm, "dm_chr_maps/dm_maps_50kb_divergence.tsv",
            quote = F, sep = "\t", col.names = T, row.names = F)

cor.test(x = lands.divergence.dm$divergence, y = lands.divergence.dm$diversity, method = "spearman")
# 0.21 p-value = 4.5e-10
cor.test(lands.divergence.dm$divergence, lands.divergence.dm$theta, method = "spearman")
# 0.20 p-value = 3e-09
cor.test(x = lands.divergence.dm$divergence, y = lands.divergence.dm$rho, method = "spearman")
# -0.003 p-value = 0.94
pcor.test(x = lands.divergence.dm$divergence, y = lands.divergence.dm$diversity,
          lands.divergence.dm$theta, method = "spearman")
# 0.07 p-value = 0.05


########################################
#
# Evolutionary (Protein) Rates --- 50 kb maps for ++ resolution
#
########################################

# loads 
dm.raw <- read.table("~/Data/iSMC/theta_paper/real_data/dm_misc/dpgp3_Dyak_bpp.all.csv", header = T, fill = T, stringsAsFactors = T)
dm.tbl <- na.omit(dm.raw)

# gets ratios (following Filipa's instructions)
dm.tbl$PiS <- dm.tbl$PiS / dm.tbl$MeanNumberSynPos
dm.tbl$PiN <- dm.tbl$PiN / (3 - dm.tbl$MeanNumberSynPos)
dm.tbl$dS <- as.numeric(dm.tbl$dS) / dm.tbl$MeanNumberSynPosDiv
dm.tbl$dN <- dm.tbl$dN / (3 - dm.tbl$MeanNumberSynPosDiv)
# cleans
dm.tbl.popgen <- as.data.frame(cbind(dm.tbl$PiN, dm.tbl$PiS, dm.tbl$dN, dm.tbl$dS, dm.tbl$GeneID))
dm.tbl.popgen <- na.omit(dm.tbl.popgen)
names(dm.tbl.popgen) <- c("PiN", "PiS", "dN", "dS", "geneID")

# for each gene, sums ratios of each codon
dm.tbl.genes <- ddply(.data = dm.tbl.popgen, .variables = "geneID", .fun = colSums, na.rm = T, .progress = "text")
# substitutes gene id and computes ratios
dm.tbl.genes$geneID <- unique(dm.tbl$GeneID)
dm.tbl.genes$dNdS <- dm.tbl.genes$dN / dm.tbl.genes$dS
dm.tbl.genes$PiNPiS <- dm.tbl.genes$PiN / dm.tbl.genes$PiS
dm.tbl.popstats <- cbind.data.frame(as.character(dm.tbl.genes$geneID), dm.tbl.genes$PiN, dm.tbl.genes$PiS, dm.tbl.genes$PiNPiS,
                                    dm.tbl.genes$dN, dm.tbl.genes$dS, dm.tbl.genes$dNdS)
names(dm.tbl.popstats) <- c("geneID", "PiN", "PiS", "PiNPiS", "dN", "dS", "dNdS")
dm.tbl.popstats$PiS <- as.numeric(dm.tbl.popstats$PiS)
dm.tbl.popstats$dS <- as.numeric(dm.tbl.popstats$dS)
dm.tbl.popstats$PiN <- as.numeric(dm.tbl.popstats$PiN)
dm.tbl.popstats$dN <- as.numeric(dm.tbl.popstats$dN)
dm.tbl.popstats$PiNPiS <- as.numeric(dm.tbl.popstats$PiNPiS)
dm.tbl.popstats$dNdS <- as.numeric(dm.tbl.popstats$dNdS)

dm.tbl.popstats.clean <- dm.tbl.popstats[which(dm.tbl.popstats$PiNPiS > 0),]
dm.tbl.popstats.clean <- dm.tbl.popstats.clean[which(dm.tbl.popstats.clean$dNdS > 0),]

dm.genes.coord <- read.table("~/Data/iSMC/theta_paper/real_data/dm_misc/Dmel_trans_nooverlaps-final.txt", header = F)
names(dm.genes.coord) <- c("chr", "start", "end", "x", "geneID", "length")
dm.genes.coord <- dm.genes.coord[,-4]

dm.evol <- merge(dm.genes.coord, dm.tbl.popstats.clean, by = "geneID")
dm.evol <- dm.evol[order(dm.evol$chr),]

# reorder maps to use GR
dm.maps.50kb <- dm.lands.50kb[,c(8,1:6)]
dm.maps.50kb$chr[which(dm.maps.50kb$chr == 2)] <- "2L"
dm.maps.50kb$chr[which(dm.maps.50kb$chr == 3)] <- "2R"
dm.maps.50kb$chr[which(dm.maps.50kb$chr == 4)] <- "3L"
dm.maps.50kb$chr[which(dm.maps.50kb$chr == 5)] <- "3R"

# grouping per gene coordinate
dm.lands.gr <- makeGRangesFromDataFrame(dm.maps.50kb)
values(dm.lands.gr) <- dm.maps.50kb[,(4:7)]
evolrate.gr <- makeGRangesFromDataFrame(dm.evol)
values(evolrate.gr) <- dm.evol[,(5:11)]

hits <- findOverlaps(query = evolrate.gr, subject = dm.lands.gr, type = "within") 

evolrate.gr.df <- as.data.frame(evolrate.gr[queryHits(hits)], row.names = NULL)
dm.lands.gr.df <- as.data.frame(dm.lands.gr[subjectHits(hits)], row.names = NULL)

dm.lands.evolrate <- cbind.data.frame(dm.lands.gr.df[,c(1:3,6:9)], evolrate.gr.df[,c(2,3,6:12)])
dm.lands.evolrate <- dm.lands.evolrate[which(dm.lands.evolrate$PiNPiS < 1),]
names(dm.lands.evolrate)[1] <- "chr"
names(dm.lands.evolrate)[2] <- "start.window"
names(dm.lands.evolrate)[3] <- "end.window"
names(dm.lands.evolrate)[8] <- "start.gene"
names(dm.lands.evolrate)[9] <- "end.gene"

# NOTE: must sort rows of dm.lands.evolrate by chr and start coords
dm.lands.evolrate[order(dm.lands.evolrate$chr),]

write.table(dm.lands.evolrate, "dm_chr_maps/dm_maps_50kb_protein_rates.tsv",
            sep = "\t", quote = F, col.names = T, row.names = F)

# linear model in coding regions
m.diversity.cds <- lm(diversity ~ (theta + rho + tmrca) * chr, data = dm.lands.evolrate)
plot(m.diversity.cds, which = 2) 
shapiro.test(resid(m.diversity.cds)) # ***
hmctest(m.diversity.cds, nsim = 3000) # NS
dwtest(m.diversity.cds) # ***

bc.diversity.cds <- boxcox(m.diversity.cds, lambda = seq(0.4, 1.1, len = 500))
l <- bc.diversity.cds$x[which.max(bc.diversity.cds$y)]
m.diversity.cds.bc <- update(m.diversity.cds, (diversity^l -1)/l~.)
plot(m.diversity.cds.bc, which = 2)
shapiro.test(resid(m.diversity.cds.bc)) # *** a bit better
hist(resid(m.diversity.cds.bc))

summary(m.diversity.cds.bc)
# Coefficients:
# Estimate Std. Error   t value Pr(>|t|)    
# (Intercept) -1.414913   0.001100 -1286.091  < 2e-16 ***
# theta        3.681421   0.057568    63.949  < 2e-16 ***
# rho          0.009793   0.012704     0.771 0.441089    
# tmrca        0.041239   0.001612    25.583  < 2e-16 ***
# chr2R       -0.003939   0.001829    -2.153 0.031724 *  
# chr3L        0.002695   0.001275     2.114 0.034967 *  
# chr3R       -0.004764   0.001299    -3.667 0.000267 ***
# theta:chr2R  0.048188   0.070104     0.687 0.492114    
# theta:chr3L  0.193858   0.064298     3.015 0.002682 ** 
# theta:chr3R  0.305470   0.068006     4.492 8.52e-06 ***
# rho:chr2R   -0.009671   0.018604    -0.520 0.603376    
# rho:chr3L   -0.001866   0.016160    -0.115 0.908100    
# rho:chr3R   -0.034404   0.020296    -1.695 0.090581 .  
# tmrca:chr2R  0.003373   0.002413     1.398 0.162658    
# tmrca:chr3L -0.005338   0.001891    -2.824 0.004910 ** 
# tmrca:chr3R  0.002920   0.002001     1.459 0.145118    

# Residual standard error: 0.0009685 on 583 degrees of freedom
# Multiple R-squared:  0.9949,	Adjusted R-squared:  0.9948

# type 2 anova
anova.diversity.cds <- Anova(m.diversity.cds.bc)
apiss <- anova.diversity.cdsc$"Sum Sq"
anova.diversity.cds$VarExp <- apiss / sum(apiss)

anova.diversity.cds
# Sum Sq  Df    F value  Pr(>F)  VarExp
# theta     0.040510   1 43184.2514 0.00000 0.91407
# rho       0.000000   1     0.0707 0.79045 0.00020
# tmrca     0.003698   1  3942.5420 0.00000 0.06050
# chr       0.000032   3    11.4592 0.00000 0.00142
# theta:chr 0.000031   3    10.8925 0.00000 0.00043
# rho:chr   0.000003   3     1.1904 0.31264 0.00023
# tmrca:chr 0.000034   3    12.1752 0.00000 0.00200
# Residuals 0.000547 583                    0.02114

# Because of auto-correlation we compute p-values for the variables using a GLS
dm.lands.evolrate$bin <- 1:nrow(dm.lands.evolrate)
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = dm.lands.evolrate, corr = corAR1(0, ~bin))

summary(g.diversity)
# Parameter estimate(s):
# Phi 
# 0.04553327 

# Coefficients:
# Value   Std.Error   t-value p-value
# (Intercept) -0.0080347 0.000145005 -55.40991  0.0000
# theta        0.9936184 0.005808070 171.07548  0.0000
# rho          0.0035568 0.002188371   1.62533  0.1046
# tmrca        0.0084526 0.000211426  39.97914  0.0000

# correlations
cor.test(dm.lands.evolrate$PiN, dm.lands.evolrate$theta, method = "spearman") 
# 0.08, p-value 0.06
pcor.test(dm.lands.evolrate$PiN, dm.lands.evolrate$theta, dm.lands.evolrate$tmrca, method = "spearman") 
# -0.0009825256 p-value 0.98
pcor.test(dm.lands.evolrate$PiS, dm.lands.evolrate$theta, dm.lands.evolrate$tmrca, method = "spearman") 
# 0.28 p-value 1.3e-12
cor.test(dm.lands.evolrate$dS, dm.lands.evolrate$PiS, method = "spearman")
# 0.65 p-value < 2.2e-16
pcor.test(dm.lands.evolrate$dS, dm.lands.evolrate$PiS, dm.lands.evolrate$theta, method = "spearman")
# 0.7 p-value 2.1e-90
cor.test(dm.lands.evolrate$dS, dm.lands.evolrate$theta, method = "spearman")
# -0.01 p-value 0.76
pcor.test(dm.lands.evolrate$dS, dm.lands.evolrate$theta, dm.lands.evolrate$tmrca, method = "spearman")
# -0.4 p-value 0.34

# checking about rec rate
cor.test(dm.lands.evolrate$PiS, dm.lands.evolrate$rho, method = "spearman") 
pcor.test(dm.lands.evolrate$PiS, dm.lands.evolrate$rho, dm.lands.evolrate$tmrca, method = "spearman") 
cor.test(dm.lands.evolrate$dS, dm.lands.evolrate$rho, method = "spearman") 

#####################
#
# B-value statistic across the Drosophila Genome
#
#####################

# TODO


