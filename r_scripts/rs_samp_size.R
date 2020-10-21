# Created: 17/06/2019
# Last modified: 12/08/2020
# Author: Gustavo Barroso


library(reshape2)
library(tidyverse)
library(scales)
library(cowplot)
library(MASS)
library(lmtest)
library(ape)
library(nlme)
library(plyr)
library(interactions)
library(rr2)

bin_sizes <- c(50000, 200e+3, 1e+6)

theme.blank <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

no.legend <- theme(legend.position = "none", legend.title = NULL)


setwd("~")
setwd("Data/iSMC/theta_paper/sim_data/rs.sim_results/")

###################################################
#
# Param estimates
#
###################################################

# parameters used in the simulation
sim_params <- read.table("raw_data/Misc/sim_params.txt")
sim_params <- sim_params[c(3:6, 8, 9),]
sim_params[4, 2] <- 1 - sim_params[4, 2]
sim_params[6, 2] <- 1 - sim_params[6, 2]
names(sim_params) <- c("param", "val")

# parameters inferred by iSMC
mle.list <- list()
for(i in 1:45) { # 45 pairs of haploids
  tmp <- readLines(paste("rep_1/hap_pair_", i, "/rs.pair_", i, "_estimates.txt", sep = ""))
  row <- c(as.numeric(strsplit(tmp[13], " ")[[1]][2]), # rho
           as.numeric(strsplit(tmp[14], " ")[[1]][2]), # theta.alpha
           as.numeric(strsplit(tmp[15], " ")[[1]][2]),  # rho.alpha
           as.numeric(strsplit(tmp[16], " ")[[1]][2]), # t_ii
           as.numeric(strsplit(tmp[17], " ")[[1]][2]), # r_ii
           as.numeric(strsplit(tmp[18], " ")[[1]][2]), # theta
           paste("pair_", i, sep = ""))
  
  mle.list[[i]] <- row
}

sim.mle <- as.data.frame(do.call(rbind, mle.list), stringsAsFactors = F)
names(sim.mle) <- c("rho", "theta.alpha", "rho.alpha", "theta.delta", "rho.delta", "theta", "sample")

# 5 "unphased" pairs of genomes
tmp <- readLines("rep_1.unphased/rep_1.unphased_estimates.txt")
row <- c(as.numeric(strsplit(tmp[13], " ")[[1]][2]), # rho
         as.numeric(strsplit(tmp[14], " ")[[1]][2]), # theta.alpha
         as.numeric(strsplit(tmp[15], " ")[[1]][2]),  # rho.alpha
         as.numeric(strsplit(tmp[16], " ")[[1]][2]), # t_ii
         as.numeric(strsplit(tmp[17], " ")[[1]][2]), # r_ii
         as.numeric(strsplit(tmp[18], " ")[[1]][2]), # theta
         "unphased")

sim.mle <- rbind(sim.mle, row)

# 45 "phased" pairs of genomes
tmp <- readLines("rep_1.joint/rep_1.joint_estimates.txt")
row <- c(as.numeric(strsplit(tmp[13], " ")[[1]][2]), # rho
         as.numeric(strsplit(tmp[14], " ")[[1]][2]), # theta.alpha
         as.numeric(strsplit(tmp[15], " ")[[1]][2]),  # rho.alpha
         as.numeric(strsplit(tmp[16], " ")[[1]][2]), # t_ii
         as.numeric(strsplit(tmp[17], " ")[[1]][2]), # r_ii
         as.numeric(strsplit(tmp[18], " ")[[1]][2]), # theta
         "phased")

sim.mle <- rbind(sim.mle, row)

# transforms to delta
sim.mle$theta.delta <- 1 - as.numeric(sim.mle$theta.delta)
sim.mle$rho.delta <- 1 - as.numeric(sim.mle$rho.delta)

sim.mle$type <- c(rep("single_diploid", 45), "5-unphased", "45-phased")

# plotting
scaleFUN <- function(x) sprintf("%.4f", x) # digits shown in y axis
no.x.axis <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())


rho.delta.df <- as.data.frame(cbind(as.numeric(sim.mle$rho.delta), sim.mle$sample, sim.mle$type), stringsAsFactors = F)
names(rho.delta.df) <- c("rho.delta", "sample", "type")
rho.delta.plot <- ggplot(rho.delta.df, aes(x = type, y = as.numeric(rho.delta), fill = type))
rho.delta.plot <- rho.delta.plot + scale_y_continuous(limits = c(0, 1.5e-5)) + no.x.axis + no.legend
rho.delta.plot <- rho.delta.plot + geom_hline(yintercept = sim_params[4,2], colour = "magenta")
rho.delta.plot <- rho.delta.plot + labs(title = NULL, x = NULL, y = "rho.delta")
rho.delta.plot <- rho.delta.plot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                binwidth = 5e-7, dotsize = 1.2, position = position_dodge(width = 0.5),
                                                alpha = 0.7, stackratio = 0.5)

rho.delta.plot

theta.delta.df <- as.data.frame(cbind(as.numeric(sim.mle$theta.delta), sim.mle$sample, sim.mle$type), stringsAsFactors = F)
names(theta.delta.df) <- c("theta.delta", "sample", "type")
theta.delta.plot <- ggplot(theta.delta.df, aes(x = type, y = as.numeric(theta.delta), fill = type))
theta.delta.plot <- theta.delta.plot + scale_y_continuous(limits = c(0, 1e-5)) + no.x.axis + no.legend
theta.delta.plot <- theta.delta.plot + geom_hline(yintercept = sim_params[6,2], colour = "magenta")
theta.delta.plot <- theta.delta.plot + labs(title = NULL, x = NULL, y = "theta.delta")
theta.delta.plot <- theta.delta.plot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                    binwidth = 5e-7, dotsize = 0.9, position = position_dodge(width = 0.5),
                                                    alpha = 0.7, stackratio = 0.5)

theta.delta.plot

theta.alpha.df <- as.data.frame(cbind(as.numeric(sim.mle$theta.alpha), sim.mle$sample, sim.mle$type), stringsAsFactors = F)
names(theta.alpha.df) <- c("theta.alpha", "sample", "type")
theta.alpha.plot <- ggplot(theta.alpha.df, aes(x = type, y = as.numeric(theta.alpha), fill = type))
theta.alpha.plot <- theta.alpha.plot + no.x.axis + scale_y_continuous(labels = scaleFUN) + no.legend
theta.alpha.plot <- theta.alpha.plot + geom_hline(yintercept = sim_params[5,2], colour = "magenta")
theta.alpha.plot <- theta.alpha.plot + labs(title = NULL, x = NULL, y = "theta.alpha")
theta.alpha.plot <- theta.alpha.plot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                    binwidth = 0.1, dotsize = 1, position = position_dodge(width = 0.5), alpha = 0.7)

theta.alpha.plot


rho.alpha.df <- as.data.frame(cbind(as.numeric(sim.mle$rho.alpha), sim.mle$sample, sim.mle$type), stringsAsFactors = F)
names(rho.alpha.df) <- c("rho.alpha", "sample", "type")
rho.alpha.plot <- ggplot(rho.alpha.df, aes(x = type, y = as.numeric(rho.alpha), fill = type))
rho.alpha.plot <- rho.alpha.plot + no.x.axis + scale_y_continuous(labels = scaleFUN) + no.legend
rho.alpha.plot <- rho.alpha.plot + geom_hline(yintercept = sim_params[3,2], colour = "magenta")
rho.alpha.plot <- rho.alpha.plot + labs(title = NULL, x = NULL, y = "rho.alpha")
rho.alpha.plot <- rho.alpha.plot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                binwidth = 0.02, dotsize = 0.45, stackratio = 0.5,
                                                position = position_dodge(width = 0.5), alpha = 0.7)

rho.alpha.plot


theta.df <- as.data.frame(cbind(as.numeric(sim.mle$theta), sim.mle$sample, sim.mle$type), stringsAsFactors = F)
names(theta.df) <- c("theta", "sample", "type")
theta.plot <- ggplot(theta.df, aes(x = type, y = as.numeric(theta), fill = type))
theta.plot <- theta.plot + ylim(0.002, 0.004) + no.x.axis + no.legend
theta.plot <- theta.plot + geom_hline(yintercept = sim_params[1,2], colour = "magenta")
theta.plot <- theta.plot + labs(title = NULL, x = "Type", y = expression(pi))
theta.plot <- theta.plot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.00002, dotsize = 4, stackratio = 0.5, alpha = 0.5)

theta.plot

rho.df <- as.data.frame(cbind(2 * as.numeric(sim.mle$rho), sim.mle$sample, sim.mle$type), stringsAsFactors = F)
names(rho.df) <- c("rho", "sample", "type")
rho.plot <- ggplot(rho.df, aes(x = type, y = as.numeric(rho), fill = type))
rho.plot <- rho.plot + scale_y_continuous(labels = scaleFUN, limits = c(0, 0.003)) + no.x.axis + no.legend
rho.plot <- rho.plot + geom_hline(yintercept = sim_params[2,2], colour = "magenta")
rho.plot <- rho.plot + labs(title = NULL, x = NULL, y = expression(rho))
rho.plot <- rho.plot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                    binwidth = 0.0002, dotsize = 0.7, position = position_dodge(width = 0.5),
                                    stackratio = 0.5, alpha = 0.7)

rho.plot


# all together now
theme_set(theme_cowplot(font_size = 12))
legend <- get_legend(rho.plot + theme(legend.position="bottom"))

p <- plot_grid(rho.plot, theta.plot,
               rho.alpha.plot, theta.alpha.plot,
               rho.delta.plot, theta.delta.plot,
               nrow = 3, ncol = 2, labels = NULL)

mle.fig <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))

cowplot::ggsave("mle.pdf", mle.fig, width = 12, height = 12)


###################################################
#
# Rho Landscapes
#
###################################################
  
# sim landscapes 50 kb
sim.rho.50k <- read.table("raw_data/rig.sims.rho.50000.bins.txt")

nbins <- 3e+7 / 5e+4
  
# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
rho.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".rho.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.rho.50k$sim, tmp$diploid_1, method = "spearman")$estimate
  rho.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.50k <- as.data.frame(cor.vec)
cor.50k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.rho.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
# 0.92
cor.test(sim.rho.50k$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.50k <- c(cor.test(sim.rho.50k$sim, unphased$sample_mean, method = "spearman")$estimate, "mean_5")
cor.50k <- rbind(cor.50k, unphased.50k)

phased <- read.table(paste("maps/rep_1.joint.rho.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
# 0.95
cor.test(sim.rho.50k$sim, phased$sample_mean, method = "spearman")$estimate
phased.50k <- c(cor.test(sim.rho.50k$sim, phased$sample_mean, method = "spearman")$estimate, "mean_45")
cor.50k <- rbind(cor.50k, phased.50k)

mean.tmrca.vec <- numeric(length = 45)
var.tmrca.vec <- numeric(length = 45)
tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
  mean.tmrca.vec[i] <- mean(tmp$diploid_1)
  var.tmrca.vec[i] <- var(tmp$diploid_1)
  tmrca.df[,i] <- tmp$diploid_1
}

# builds "consensus" rec. map by cherry-picking local TMRCA values
pair.ids.high.tmrca <- apply(tmrca.df, 1, which.max)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.max <- pair.ids.high.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.max]
}
# 0.89
cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate, "highest.tmrca")
cor.50k <- rbind(cor.50k, high.tmrca)

# taking the mean rho among top 5 TMRCA values
top.5.tmrca <- list()
for(i in 1:nrow(tmrca.df)) {
  tmp <- sort(tmrca.df[i,])
  top.5.tmrca[[i]] <- tmp[41:45]
}
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  top.5.rho <- numeric(length = 5)
  max.pairs <- top.5.tmrca[[i]]
  tmp <- rho.df[i,]
  for(j in 1:5) {
    top.5.rho[j] <- as.numeric(tmp[names(max.pairs)[j]])
  }
  cherry.pick.rec.map[i] <- mean(top.5.rho)
}
# 0.93
cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate, "top5.tmrca")
cor.50k <- rbind(cor.50k, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.min]
}
# 0.57
cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate 
low.tmrca <- c(cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate, "lowest.tmrca")
cor.50k <- rbind(cor.50k, low.tmrca)

# preparing to plot
cor.50k$cor.vec <- as.numeric(cor.50k$cor.vec)
cor.50k$type <- factor(cor.50k$type)
cor.50k$type <- factor(cor.50k$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

cor.dotplot <- ggplot(data = cor.50k, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                          binwidth = 0.005, dotsize = 4, position = position_dodge(width=0.1), alpha = 0.8)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1) + theme.blank
cor.dotplot





# 200 kb

sim.rho.200k <- read.table("raw_data/rig.sims.rho.2e+05.bins.txt")

nbins <- 3e+7 / 2e+5

# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
rho.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".rho.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.rho.200k$sim, tmp$diploid_1, method = "spearman")$estimate
  rho.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.200k <- as.data.frame(cor.vec)
cor.200k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.rho.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
# 0.93
cor.test(sim.rho.200k$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.200k <- c(cor.test(sim.rho.200k$sim, unphased$sample_mean, method = "spearman")$estimate, "mean_5")
cor.200k <- rbind(cor.200k, unphased.200k)

phased <- read.table(paste("maps/rep_1.joint.rho.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
# 0.95
cor.test(sim.rho.200k$sim, phased$sample_mean, method = "spearman")$estimate
phased.200k <- c(cor.test(sim.rho.200k$sim, unphased$sample_mean, method = "spearman")$estimate, "mean_45")
cor.200k <- rbind(cor.200k, phased.200k)


tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
  tmrca.df[,i] <- tmp$diploid_1
}

# builds "consensus" rec. map by cherry-picking based on local TMRCA values
pair.ids.high.tmrca <- apply(tmrca.df, 1, which.max)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.max <- pair.ids.high.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.max]
}
# 0.87
cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate, "highest.tmrca")
cor.200k <- rbind(cor.200k, high.tmrca)

# taking the mean rho among top 5 TMRCA values
top.5.tmrca <- list()
for(i in 1:nrow(tmrca.df)) {
  tmp <- sort(tmrca.df[i,])
  top.5.tmrca[[i]] <- tmp[41:45]
}
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  top.5.rho <- numeric(length = 5)
  max.pairs <- top.5.tmrca[[i]]
  tmp <- rho.df[i,]
  for(j in 1:5) {
    top.5.rho[j] <- as.numeric(tmp[names(max.pairs)[j]])
  }
  cherry.pick.rec.map[i] <- mean(top.5.rho)
}
# 0.91
cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate, "top5.tmrca")
cor.200k <- rbind(cor.200k, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.min]
}
# 0.60
cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate, "lowest.tmrca")
cor.200k <- rbind(cor.200k, low.tmrca)

# preparing to plot
cor.200k$cor.vec <- as.numeric(cor.200k$cor.vec)
cor.200k$type <- factor(cor.200k$type)
cor.200k$type <- factor(cor.200k$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

cor.dotplot <- ggplot(data = cor.200k, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.5,
                                        binwidth = 0.005, dotsize = 4, position = position_dodge(width=0.5), alpha = 0.8)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman Correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1) + theme.blank
cor.dotplot





# 1 Mb

sim.rho.1M <- read.table("raw_data/rig.sims.rho.1e+06.bins.txt")

nbins <- 3e+7 / 1e+6

# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
rho.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".rho.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.rho.1M$sim, tmp$diploid_1, method = "spearman")$estimate
  rho.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.1M <- as.data.frame(cor.vec)
cor.1M$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.rho.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
# 0.90
cor.test(sim.rho.1M$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.1M <- c(cor.test(sim.rho.1M$sim, unphased$sample_mean, method = "spearman")$estimate, "mean_5")
cor.1M <- rbind(cor.1M, unphased.1M)

phased <- read.table(paste("maps/rep_1.joint.rho.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
# 0.92
cor.test(sim.rho.1M$sim, phased$sample_mean, method = "spearman")$estimate
phased.1M <- c(cor.test(sim.rho.1M$sim, phased$sample_mean, method = "spearman")$estimate, "mean_45")
cor.1M <- rbind(cor.1M, phased.1M)

tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
  tmrca.df[,i] <- tmp$diploid_1
}

# builds "consensus" rec. map by cherry-picking based on local TMRCA values
pair.ids.high.tmrca <- apply(tmrca.df, 1, which.max)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.max <- pair.ids.high.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.max]
}
# 0.89
cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate, "highest.tmrca")
cor.1M <- rbind(cor.1M, high.tmrca)

# taking the mean rho among top 5 TMRCA values
top.5.tmrca <- list()
for(i in 1:nrow(tmrca.df)) {
  tmp <- sort(tmrca.df[i,])
  top.5.tmrca[[i]] <- tmp[41:45]
}
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  top.5.rho <- numeric(length = 5)
  max.pairs <- top.5.tmrca[[i]]
  tmp <- rho.df[i,]
  for(j in 1:5) {
    top.5.rho[j] <- as.numeric(tmp[names(max.pairs)[j]])
  }
  cherry.pick.rec.map[i] <- mean(top.5.rho)
}
# 0.89
cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate, "top5.tmrca")
cor.1M <- rbind(cor.1M, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.min]
}
# 0.79
cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate, "lowest.tmrca")
cor.1M <- rbind(cor.1M, low.tmrca)

# preparing to plot
cor.1M$cor.vec <- as.numeric(cor.1M$cor.vec)
cor.1M$type <- factor(cor.1M$type)
cor.1M$type <- factor(cor.1M$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

cor.dotplot <- ggplot(data = cor.1M, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.005, dotsize = 4, position = position_dodge(width=0.5), alpha = 0.75)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1) + theme.blank
cor.dotplot

# all together now
cor.df <- rbind(cor.50k, cor.200k, cor.1M)
cor.df$scale <- c(rep("50kb", 50), rep("200kb", 50), rep("1Mb", 50))

cor.dotplot <- ggplot(data = cor.df, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.0005, dotsize = 40, position = position_dodge(width=0.7), alpha = 0.75) + ylim(0, 1)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "cor Rec. Map") + theme.blank
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~scale) 
ggplot2::ggsave("cor.rho.pdf", cor.dotplot, device = "pdf", width = 12, height = 9)



###################################################
#
# Theta Landscapes
#
###################################################


# 50 kb
sim.theta.50k <- read.table("raw_data/rig.sims.theta.50000.bins.txt")

nbins <- 3e+7 / 5e+4

# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
theta.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".theta.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.theta.50k$sim, tmp$diploid_1, method = "spearman")$estimate
  theta.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.50k <- as.data.frame(cor.vec)
cor.50k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.theta.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
# 0.90
cor.test(sim.theta.50k$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.50k <- c(cor.test(sim.theta.50k$sim, unphased$sample_mean, method = "spearman")$estimate, "mean_5")
cor.50k <- rbind(cor.50k, unphased.50k)

phased <- read.table(paste("maps/rep_1.joint.theta.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
# 0.90
cor.test(sim.theta.50k$sim, phased$sample_mean, method = "spearman")$estimate
phased.50k <- c(cor.test(sim.theta.50k$sim, phased$sample_mean, method = "spearman")$estimate, "mean_45")
cor.50k <- rbind(cor.50k, phased.50k)

tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
  tmrca.df[,i] <- tmp$diploid_1
}

# builds "consensus" mut. map by cherry-picking local TMRCA values
pair.ids.high.tmrca <- apply(tmrca.df, 1, which.max)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.max <- pair.ids.high.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.max]
}
# 0.88
cor.test(sim.theta.50k$sim, cherry.pick.mut.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.theta.50k$sim, cherry.pick.mut.map, method = "spearman")$estimate, "highest.tmrca")
cor.50k <- rbind(cor.50k, high.tmrca)

# taking the mean rho among top 5 TMRCA values
top.5.tmrca <- list()
for(i in 1:nrow(tmrca.df)) {
  tmp <- sort(tmrca.df[i,])
  top.5.tmrca[[i]] <- tmp[41:45]
}
cherry.pick.rec.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  top.5.rho <- numeric(length = 5)
  max.pairs <- top.5.tmrca[[i]]
  tmp <- theta.df[i,]
  for(j in 1:5) {
    top.5.rho[j] <- as.numeric(tmp[names(max.pairs)[j]])
  }
  cherry.pick.rec.map[i] <- mean(top.5.rho)
}
# 0.89
cor.test(sim.theta.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.theta.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate, "top5.tmrca")
cor.50k <- rbind(cor.50k, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.min]
}
# 0.60
cor.test(sim.theta.50k$sim, cherry.pick.mut.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.theta.50k$sim, cherry.pick.mut.map, method = "spearman")$estimate, "lowest.tmrca")
cor.50k <- rbind(cor.50k, low.tmrca)

# preparing to plot
cor.50k$cor.vec <- as.numeric(cor.50k$cor.vec)
cor.50k$type <- factor(cor.50k$type)
cor.50k$type <- factor(cor.50k$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

cor.dotplot <- ggplot(data = cor.50k, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.02, dotsize = 1, position = position_dodge(width=0.5), alpha = 0.75)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1) + theme.blank
cor.dotplot





# 200 kb

sim.theta.200k <- read.table("raw_data/rig.sims.theta.2e+05.bins.txt")

nbins <- 3e+7 / 2e+5

# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
theta.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".theta.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.theta.200k$sim, tmp$diploid_1, method = "spearman")$estimate
  theta.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.200k <- as.data.frame(cor.vec)
cor.200k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.theta.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
# 0.92
cor.test(sim.theta.200k$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.200k <- c(cor.test(sim.theta.200k$sim, unphased$sample_mean, method = "spearman")$estimate, "mean_5")
cor.200k <- rbind(cor.200k, unphased.200k)

phased <- read.table(paste("maps/rep_1.joint.theta.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
# 0.93
cor.test(sim.theta.200k$sim, phased$sample_mean, method = "spearman")$estimate 
phased.200k <- c(cor.test(sim.theta.200k$sim, phased$sample_mean, method = "spearman")$estimate, "mean_45")
cor.200k <- rbind(cor.200k, phased.200k)

tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
  tmrca.df[,i] <- tmp$diploid_1
}

# builds "consensus" mut. map by cherry-picking local TMRCA values
pair.ids.high.tmrca <- apply(tmrca.df, 1, which.max)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.max <- pair.ids.high.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.max]
}
# 0.88
cor.test(sim.theta.200k$sim, cherry.pick.mut.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.theta.200k$sim, cherry.pick.mut.map, method = "spearman")$estimate, "highest.tmrca")
cor.200k <- rbind(cor.200k, high.tmrca)

# taking the mean rho among top 5 TMRCA values
top.5.tmrca <- list()
for(i in 1:nrow(tmrca.df)) {
  tmp <- sort(tmrca.df[i,])
  top.5.tmrca[[i]] <- tmp[41:45]
}
cherry.pick.rec.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  top.5.rho <- numeric(length = 5)
  max.pairs <- top.5.tmrca[[i]]
  tmp <- theta.df[i,]
  for(j in 1:5) {
    top.5.rho[j] <- as.numeric(tmp[names(max.pairs)[j]])
  }
  cherry.pick.rec.map[i] <- mean(top.5.rho)
}
# 0.93
cor.test(sim.theta.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.theta.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate, "top5.tmrca")
cor.200k <- rbind(cor.200k, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.min]
}
# 0.58
cor.test(sim.theta.200k$sim, cherry.pick.mut.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.theta.200k$sim, cherry.pick.mut.map, method = "spearman")$estimate, "lowest.tmrca")
cor.200k <- rbind(cor.200k, low.tmrca)

# preparing to plot
cor.200k$cor.vec <- as.numeric(cor.200k$cor.vec)
cor.200k$type <- factor(cor.200k$type)
cor.200k$type <- factor(cor.200k$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

cor.dotplot <- ggplot(data = cor.200k, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.005, dotsize = 4, position = position_dodge(width=0.5), alpha = 0.75)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1) + theme.blank
cor.dotplot




# 1 Mb

sim.theta.1M <- read.table("raw_data/rig.sims.theta.1e+06.bins.txt")

nbins <- 3e+7 / 1e+6

# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
theta.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".theta.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.theta.1M$sim, tmp$diploid_1, method = "spearman")$estimate
  theta.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.1M <- as.data.frame(cor.vec)
cor.1M$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.theta.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
# 0.93
cor.test(sim.theta.1M$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.1M <- c(cor.test(sim.theta.1M$sim, unphased$sample_mean, method = "spearman")$estimate, "mean_5")
cor.1M <- rbind(cor.1M, unphased.1M)

phased <- read.table(paste("maps/rep_1.joint.theta.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
# 0.94
cor.test(sim.theta.1M$sim, phased$sample_mean, method = "spearman")$estimate
phased.1M <- c(cor.test(sim.theta.1M$sim, phased$sample_mean, method = "spearman")$estimate, "mean_45")
cor.1M <- rbind(cor.1M, phased.1M)

tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
  tmrca.df[,i] <- tmp$diploid_1
}

# builds "consensus" mut. map by cherry-picking local TMRCA values
pair.ids.high.tmrca <- apply(tmrca.df, 1, which.max)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.max <- pair.ids.high.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.max]
}
# 0.91
cor.test(sim.theta.1M$sim, cherry.pick.mut.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.theta.1M$sim, cherry.pick.mut.map, method = "spearman")$estimate, "highest.tmrca")
cor.1M <- rbind(cor.1M, high.tmrca)

# taking the mean rho among top 5 TMRCA values
top.5.tmrca <- list()
for(i in 1:nrow(tmrca.df)) {
  tmp <- sort(tmrca.df[i,])
  top.5.tmrca[[i]] <- tmp[41:45]
}
cherry.pick.rec.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  top.5.rho <- numeric(length = 5)
  max.pairs <- top.5.tmrca[[i]]
  tmp <- theta.df[i,]
  for(j in 1:5) {
    top.5.rho[j] <- as.numeric(tmp[names(max.pairs)[j]])
  }
  cherry.pick.rec.map[i] <- mean(top.5.rho)
}
# 0.93
cor.test(sim.theta.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.theta.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate, "top5.tmrca")
cor.1M <- rbind(cor.1M, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.min]
}
# 0.74
cor.test(sim.theta.1M$sim, cherry.pick.mut.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.theta.1M$sim, cherry.pick.mut.map, method = "spearman")$estimate, "lowest.tmrca")
cor.1M <- rbind(cor.1M, low.tmrca)

# preparing to plot
cor.1M$cor.vec <- as.numeric(cor.1M$cor.vec)
cor.1M$type <- factor(cor.1M$type)
cor.1M$type <- factor(cor.1M$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

cor.dotplot <- ggplot(data = cor.1M, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.02, dotsize = 1, position = position_dodge(width=0.5), alpha = 0.75)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1) + theme.blank
cor.dotplot

# all together now
cor.df <- rbind(cor.50k, cor.200k, cor.1M)
cor.df$scale <- c(rep("50kb", 50), rep("200kb", 50), rep("1Mb", 50))
cor.dotplot <- ggplot(data = cor.df, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.001, dotsize = 15, position = position_dodge(width=0.3), alpha = 0.75) + ylim(0, 1)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "cor Mut. Map") + theme.blank
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~scale) 
ggplot2::ggsave("cor.theta.pdf", cor.dotplot, device = "pdf", width = 12, height = 9)



###################################################
#
# TMRCA Landscapes
#
###################################################
  
# gets nucleotide spans of each tree
arg <- readLines("raw_data/rep_1/rep_1.ARG")
tmp <- list()
for(i in 1:length(arg)) {
  tmp[i] <- data.matrix(strsplit(arg[i], '"', fixed = F))
}
arg <- list()
for(i in 1:length(tmp)) {
  arg[i] <- tmp[[i]][4]
}
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

sim.ARG <- read.tree("raw_data/rep_1/rep_1.ARG")

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


nbins <- 3e+7 / 5e+4

# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
TMRCA.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.tmrca.50k$tmrca, tmp$diploid_1, method = "spearman")$estimate
  TMRCA.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.50k <- as.data.frame(cor.vec)
cor.50k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.TMRCA.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
unphased$avg <- apply(unphased[4:ncol(unphased)], 1, mean)
cor.test(sim.tmrca.50k$tmrca, unphased$avg, method = "spearman")$estimate
unphased.50k <- c(cor.test(sim.tmrca.50k$tmrca, unphased$avg, method = "spearman")$estimate, "mean_het_5")
cor.50k <- rbind(cor.50k, unphased.50k)

phased <- read.table(paste("maps/rep_1.joint.TMRCA.block.1.50kb.0-29999999.bedgraph", sep = ""), header = T)
phased$avg <- apply(phased[4:ncol(phased)], 1, mean)
# 0.67
cor.test(sim.tmrca.50k$tmrca, phased$avg, method = "spearman")$estimate
phased.50k <- c(cor.test(sim.tmrca.50k$tmrca, phased$avg, method = "spearman")$estimate, "mean_het_45")
cor.50k <- rbind(cor.50k, phased.50k)

# loads TMRCAs from 30x1x1 configuration to compare performance
tmrca.30x1x1 <- read.table("30x1x1/rep_1.unphased.30x1x1/rep_1.unphased.TMRCA.50kb.0-29999999.bedgraph", skip = 1)
tmrca.30x1x1$avg <- apply(tmrca.30x1x1[4:ncol(tmrca.30x1x1)], 1, mean)
cor.test(sim.tmrca.50k$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate
homo.50k <- c(cor.test(sim.tmrca.50k$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate, "mean_homo_5")
cor.50k <- rbind(cor.50k, homo.50k)

tmrca.30x1x1 <- read.table("30x1x1/rep_1.joint.30x1x1/rep_1.joint.TMRCA.50kb.0-29999999.bedgraph", skip = 1)
tmrca.30x1x1$avg <- apply(tmrca.30x1x1[4:ncol(tmrca.30x1x1)], 1, mean)
cor.test(sim.tmrca.50k$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate
homo.50k <- c(cor.test(sim.tmrca.50k$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate, "mean_homo_45")
cor.50k <- rbind(cor.50k, homo.50k)


# preparing to plot
cor.50k$cor.vec <- as.numeric(cor.50k$cor.vec)
cor.50k$type <- factor(cor.50k$type)
cor.50k$type <- factor(cor.50k$type, levels = c("single_diploid", "mean_homo_5", "mean_het_5", "mean_homo_45", "mean_het_45"))

cor.dotplot <- ggplot(data = cor.50k, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.02, dotsize = 1.2, position = position_dodge(width=0.5), alpha = 0.7)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)





# 200kb
tmrca.single.nuc$bin <- ceiling(1:nrow(tmrca.single.nuc) / 200e+3)
sim.tmrca.200k <- ddply(.data = tmrca.single.nuc[-which(names(tmrca.single.nuc) == "pos")],
                        .variables = "bin", .fun = colMeans, .progress = "text")


nbins <- 3e+7 / 2e+5

# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
TMRCA.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.tmrca.200k$tmrca, tmp$diploid_1, method = "spearman")$estimate
  TMRCA.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.200k <- as.data.frame(cor.vec)
cor.200k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.TMRCA.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
unphased$avg <- apply(unphased[4:ncol(unphased)], 1, mean)
cor.test(sim.tmrca.200k$tmrca, unphased$avg, method = "spearman")$estimate
unphased.200k <- c(cor.test(sim.tmrca.200k$tmrca, unphased$avg, method = "spearman")$estimate, "mean_het_5")
cor.200k <- rbind(cor.200k, unphased.200k)

phased <- read.table(paste("maps/rep_1.joint.TMRCA.block.1.200kb.0-29999999.bedgraph", sep = ""), header = T)
phased$avg <- apply(phased[4:ncol(phased)], 1, mean)
# 0.54
cor.test(sim.tmrca.200k$tmrca, phased$avg, method = "spearman")$estimate
phased.200k <- c(cor.test(sim.tmrca.200k$tmrca, phased$avg, method = "spearman")$estimate, "mean_het_45")
cor.200k <- rbind(cor.200k, phased.200k)

# loads TMRCAs from 30x1x1 configuration to compare performance
tmrca.30x1x1 <- read.table("30x1x1/rep_1.unphased.30x1x1/rep_1.unphased.TMRCA.200kb.0-29999999.bedgraph", skip = 1)
tmrca.30x1x1$avg <- apply(tmrca.30x1x1[4:ncol(tmrca.30x1x1)], 1, mean)
cor.test(sim.tmrca.200k$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate
homo.200k <- c(cor.test(sim.tmrca.200k$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate, "mean_homo_5")
cor.200k <- rbind(cor.200k, homo.200k)

tmrca.30x1x1 <- read.table("30x1x1/rep_1.joint.30x1x1/rep_1.joint.TMRCA.200kb.0-29999999.bedgraph", skip = 1)
tmrca.30x1x1$avg <- apply(tmrca.30x1x1[4:ncol(tmrca.30x1x1)], 1, mean)
cor.test(sim.tmrca.200k$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate
homo.200k <- c(cor.test(sim.tmrca.200k$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate, "mean_homo_45")
cor.200k <- rbind(cor.200k, homo.200k)


# preparing to plot
cor.200k$cor.vec <- as.numeric(cor.200k$cor.vec)
cor.200k$type <- factor(cor.200k$type)
cor.200k$type <- factor(cor.200k$type, levels = c("single_diploid", "mean_homo_5", "mean_het_5", "mean_homo_45", "mean_het_45"))

cor.dotplot <- ggplot(data = cor.200k, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.02, dotsize = 1.2, position = position_dodge(width=0.5), alpha = 0.7)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)







# 1Mb
tmrca.single.nuc$bin <- ceiling(1:nrow(tmrca.single.nuc) / 1e+6)
sim.tmrca.1M <- ddply(.data = tmrca.single.nuc[-which(names(tmrca.single.nuc) == "pos")],
                      .variables = "bin", .fun = colMeans, .progress = "text")

nbins <- 3e+7 / 1e+6

# 45 single pairs of genomes
cor.vec <- numeric(length = 45)
TMRCA.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
  cor.vec[i] <- cor.test(sim.tmrca.1M$tmrca, tmp$diploid_1, method = "spearman")$estimate
  TMRCA.df[,i] <- tmp$diploid_1
}
summary(cor.vec)

# to organise
cor.1M <- as.data.frame(cor.vec)
cor.1M$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.TMRCA.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
unphased$avg <- apply(unphased[4:ncol(unphased)], 1, mean)
cor.test(sim.tmrca.1M$tmrca, unphased$avg, method = "spearman")$estimate
unphased.1M <- c(cor.test(sim.tmrca.1M$tmrca, unphased$avg, method = "spearman")$estimate, "mean_het_5")
cor.1M <- rbind(cor.1M, unphased.1M)

phased <- read.table(paste("maps/rep_1.joint.TMRCA.block.1.1Mb.0-29999999.bedgraph", sep = ""), header = T)
phased$avg <- apply(phased[4:ncol(phased)], 1, mean)
cor.test(sim.tmrca.1M$tmrca, phased$avg, method = "spearman")$estimate
phased.1M <- c(cor.test(sim.tmrca.1M$tmrca, phased$avg, method = "spearman")$estimate, "mean_het_45")
cor.1M <- rbind(cor.1M, phased.1M)


# loads TMRCAs from 30x1x1 configuration to compare performance
tmrca.30x1x1 <- read.table("30x1x1/rep_1.unphased.30x1x1/rep_1.unphased.TMRCA.1Mb.0-29999999.bedgraph", skip = 1)
tmrca.30x1x1$avg <- apply(tmrca.30x1x1[4:ncol(tmrca.30x1x1)], 1, mean)
cor.test(sim.tmrca.1M$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate
homo.1M <- c(cor.test(sim.tmrca.1M$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate, "mean_homo_5")
cor.1M <- rbind(cor.1M, homo.1M)

tmrca.30x1x1 <- read.table("30x1x1/rep_1.joint.30x1x1/rep_1.joint.TMRCA.1Mb.0-29999999.bedgraph", skip = 1)
tmrca.30x1x1$avg <- apply(tmrca.30x1x1[4:ncol(tmrca.30x1x1)], 1, mean)
cor.test(sim.tmrca.1M$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate
homo.1M <- c(cor.test(sim.tmrca.1M$tmrca, tmrca.30x1x1$avg, method = "spearman")$estimate, "mean_homo_45")
cor.1M <- rbind(cor.1M, homo.1M)


# preparing to plot
cor.1M$cor.vec <- as.numeric(cor.1M$cor.vec)
cor.1M$type <- factor(cor.1M$type)
cor.1M$type <- factor(cor.1M$type, levels = c("single_diploid", "mean_homo_5", "mean_het_5", "mean_homo_45", "mean_het_45"))

cor.dotplot <- ggplot(data = cor.1M, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.02, dotsize = 1.2, position = position_dodge(width=0.5), alpha = 0.7)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation")
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)


# all together now
cor.df <- rbind(cor.50k, cor.200k, cor.1M)
cor.df[which(cor.df$cor.vec < 0), 1] <- 0
cor.df$scale <- c(rep("50kb", 49), rep("200kb", 49), rep("1Mb", 49))
cor.dotplot <- ggplot(data = cor.df, aes(y = cor.vec, x = type, fill = type))
cor.dotplot <- cor.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.003, dotsize = 5, position = position_dodge(width=0.3), alpha = 0.75) + ylim(0, 1)
cor.dotplot <- cor.dotplot + labs(title = NULL, x = "Map Source", y = "Spearman correlation") + theme.blank
cor.dotplot <- cor.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~scale) 
ggplot2::ggsave("cor.TMRCA.pdf", cor.dotplot, device = "pdf", width = 12, height = 9)



###################################################
#
# diversity ~ rho + theta + tmrca --- 50 kb scale (PHASED)
#
###################################################




# R_2 table for plotting at the end (for all PHASED bin sizes as well as sim.lands)
cor.tab <- as.data.frame(matrix(ncol = 5, nrow = 6))
colnames(cor.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Bin_size(kb)")




# sim landscapes 50 kb
sim.rho.50k <- read.table("raw_data/rig.sims.rho.50000.bins.txt")
sim.theta.50k <- read.table("raw_data/rig.sims.theta.50000.bins.txt")

# loading computed diversity landscape
joint.diversity.50k <- read.table("maps/rep_1.joint.diversity.block.1.50kb.0-29999999.bedgraph", header = T)
joint.diversity.50k$avg <- apply(joint.diversity.50k[4:ncol(joint.diversity.50k)], 1, mean)

# loading inferred landscapes
joint.tmrca.50k <- read.table("maps/rep_1.joint.TMRCA.block.1.50kb.0-29999999.bedgraph", header = T)
joint.tmrca.50k$avg <- apply(joint.tmrca.50k[4:ncol(joint.tmrca.50k)], 1, mean)
joint.rho.50k <- read.table("maps/rep_1.joint.rho.block.1.50kb.0-29999999.bedgraph", header = T)
joint.theta.50k <- read.table("maps/rep_1.joint.theta.block.1.50kb.0-29999999.bedgraph", header = T)

inf.lands <- as.data.frame(cbind(joint.diversity.50k$avg, joint.theta.50k$sample_mean, joint.rho.50k$sample_mean, joint.tmrca.50k$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 50, y = value, colour = "#F8766D")) + theme.blank
diversity.map <- diversity.map + geom_line(data = molten.diversity)
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, 2 * inf.lands$rho, sim.rho.50k$sim)) # adds simulated rec. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 50, y = value, colour = variable)) + theme.blank
rho.map <- rho.map + geom_line(data = molten.rho) + scale_fill_manual(values = c("#fc8d59", "#99d594"))
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.009, label = "Cor = 0.95")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.50k$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 50, y = value, colour = variable)) + theme.blank
theta.map <- theta.map + geom_line(data = molten.theta) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.0095, label = "Cor = 0.90")

tmrca.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$tmrca, sim.tmrca.50k$tmrca))
names(tmrca.plot) <- c("bin", "inf", "sim")
molten.tmrca <- melt(tmrca.plot, id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 50, y = value, colour = variable)) + theme.blank
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F) + no.legend
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(text = element_text(size = 20), title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
tmrca.map <- tmrca.map + annotate("text", x = 15000, y = 2.2, label = "Cor = 0.67")

# all together now
theme_set(theme_cowplot(font_size = 12))
legend <- get_legend(rho.map + theme(legend.position="bottom"))
p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
lands.plot.50k <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))

cowplot::save_plot("landscapes.50kb.phased.pdf", lands.plot.50k, base_width = 9, base_height = 15)

# Linear model

# centering
inf.lands$thetaC <- inf.lands$theta - mean(inf.lands$theta)
inf.lands$tmrcaC <- inf.lands$tmrca - mean(inf.lands$tmrca)
inf.lands$rhoC <- inf.lands$rho - mean(inf.lands$rho)


plot(diversity~theta, data = inf.lands)
plot(diversity~tmrca, data = inf.lands)
plot(diversity~rho, data = inf.lands)  

plot(theta~tmrca, data = inf.lands)
plot(tmrca~rho, data = inf.lands)
plot(theta~rho, data = inf.lands)

m.diversity <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands)

plot(y=resid(m.diversity, type = "pearson"), x=fitted(m.diversity))
dwtest(m.diversity)
hmctest(m.diversity)
hist(resid(m.diversity))

summary(m.diversity)
# (Intercept)   3.254e-03  1.050e-05  309.81  < 2e-16 ***
# thetaC        9.275e-01  6.357e-03  145.90  < 2e-16 ***
# rhoC          9.024e-02  1.876e-02    4.81 1.91e-06 ***
# tmrcaC        4.077e-03  4.379e-05   93.09  < 2e-16 ***
# thetaC:tmrcaC 1.218e+00  2.514e-02   48.46  < 2e-16 ***

# type 2 anova
anova.diversity <- Anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity

cor.tab[1,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                  anova.diversity$VarExp[1] * 100, 0, anova.diversity$VarExp[3] * 100, 50)

g.diversity1 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, 
                    data = inf.lands, weights = varPower(0, ~tmrcaC), corr = corAR1(0, ~bin), method = "ML")

g.diversity2 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, 
                    data = inf.lands, corr = corAR1(0, ~bin), method = "ML")

AIC(g.diversity1, g.diversity2)

summary(g.diversity2)

vif(g.diversity2)
# thetaC          rhoC        tmrcaC thetaC:tmrcaC 
# 1.107302      1.022703      1.103215      1.040263

###################################################
#
# diversity ~ rho + theta + tmrca --- 200 kb scale (PHASED)
#
###################################################

# sim landscapes 200 kb
sim.rho.200k <- read.table("raw_data/rig.sims.rho.2e+05.bins.txt")
sim.theta.200k <- read.table("raw_data/rig.sims.theta.2e+05.bins.txt")

# loading computed diversity landscape
joint.diversity.200k <- read.table("maps/rep_1.joint.diversity.block.1.200kb.0-29999999.bedgraph", header = T)
joint.diversity.200k$avg <- apply(joint.diversity.200k[4:ncol(joint.diversity.200k)], 1, mean)

# loading inferred landscapes
joint.tmrca.200k <- read.table("maps/rep_1.joint.TMRCA.block.1.200kb.0-29999999.bedgraph", header = T)
joint.tmrca.200k$avg <- apply(joint.tmrca.200k[4:ncol(joint.tmrca.200k)], 1, mean)
joint.rho.200k <- read.table("maps/rep_1.joint.rho.block.1.200kb.0-29999999.bedgraph", header = T)
joint.theta.200k <- read.table("maps/rep_1.joint.theta.block.1.200kb.0-29999999.bedgraph", header = T)

inf.lands <- as.data.frame(cbind(joint.diversity.200k$avg, joint.theta.200k$sample_mean, joint.rho.200k$sample_mean, joint.tmrca.200k$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

# does same analyses using the TRUE, SIMULATED landscapes
sim.lands <- as.data.frame(cbind(joint.diversity.200k$avg, sim.theta.200k$sim, sim.rho.200k$sim, sim.tmrca.200k$tmrca))
names(sim.lands) <- c("diversity", "theta", "rho", "tmrca")

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 200, y = value, colour = "#F8766D")) + theme.blank
diversity.map <- diversity.map + geom_line(data = molten.diversity) 
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, 2 * inf.lands$rho, sim.rho.200k$sim)) # adds simulated rec. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 200, y = value, colour = variable)) + theme.blank
rho.map <- rho.map + geom_line(data = molten.rho) 
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.007, label = "Cor = 0.95")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.200k$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 200, y = value, colour = variable)) + theme.blank
theta.map <- theta.map + geom_line(data = molten.theta)
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.0085, label = "Cor = 0.93")

tmrca.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$tmrca, sim.tmrca.200k$tmrca))
names(tmrca.plot) <- c("bin", "inf", "sim")
molten.tmrca <- melt(tmrca.plot, id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 200, y = value, colour = variable)) + theme.blank
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F) + no.legend
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
tmrca.map <- tmrca.map + annotate("text", x = 15000, y = 1.65, label = "Cor = 0.56")

# all together now
theme_set(theme_cowplot(font_size = 12))
legend <- get_legend(rho.map + theme(legend.position="bottom"))
p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
lands.plot.200k <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::save_plot("landscapes.200kb.phased.pdf", lands.plot.200k, base_width = 9, base_height = 15)



# Linear model 
inf.lands$thetaC <- inf.lands$theta - mean(inf.lands$theta)
inf.lands$tmrcaC <- inf.lands$tmrca - mean(inf.lands$tmrca)
inf.lands$rhoC <- inf.lands$rho - mean(inf.lands$rho)


plot(diversity~theta, data = inf.lands)
plot(diversity~tmrca, data = inf.lands)
plot(diversity~rho, data = inf.lands)  

plot(theta~tmrca, data = inf.lands)
plot(tmrca~rho, data = inf.lands)
plot(theta~rho, data = inf.lands)

m.diversity <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands)

plot(y=resid(m.diversity, type = "pearson"), x=fitted(m.diversity))
dwtest(m.diversity)
hmctest(m.diversity)
hist(resid(m.diversity))

summary(m.diversity)
# (Intercept)   3.310e-03  1.760e-05 188.038   <2e-16 ***
# thetaC        9.402e-01  1.150e-02  81.739   <2e-16 ***
# rhoC          6.181e-02  3.848e-02   1.606     0.11    
# tmrcaC        3.959e-03  9.998e-05  39.596   <2e-16 ***
# thetaC:tmrcaC 1.151e+00  6.188e-02  18.597   <2e-16 ***

# type 2 anova
anova.diversity <- Anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# thetaC        3.1068e-04   1 7272.7597 0.00000 0.77112
# rhoC          1.1000e-07   1    2.5807 0.11035 0.00027
# tmrcaC        7.1136e-05   1 1665.2548 0.00000 0.17656
# thetaC:tmrcaC 1.4774e-05   1  345.8588 0.00000 0.03667
# Residuals     6.1940e-06 145                   0.01537

cor.tab[2,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                 anova.diversity$VarExp[1] * 100, 0, anova.diversity$VarExp[3] * 100, 200)

g.diversity <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC,
                   data = inf.lands, weights = varPower(0, ~tmrcaC), corr = corAR1(0, ~bin))

summary(g.diversity)
# (Intercept)   0.0032884 0.00002155 152.61775  0.0000
# thetaC        0.9288933 0.01244186  74.65874  0.0000
# rhoC          0.0490388 0.03528678   1.38972  0.1667
# tmrcaC        0.0039943 0.00011366  35.14187  0.0000
# thetaC:tmrcaC 1.0940166 0.06580156  16.62600  0.0000

vif(g.diversity)
# thetaC          rhoC        tmrcaC thetaC:tmrcaC
# 1.098983      1.054555      1.085387      1.049645

###################################################
#
# diversity ~ rho + theta + tmrca --- 1Mb scale (PHASED)
#
###################################################

# sim landscapes
sim.rho.1M <- read.table("raw_data/rig.sims.rho.1e+06.bins.txt")
sim.theta.1M <- read.table("raw_data/rig.sims.theta.1e+06.bins.txt")

# loading computed diversity landscape
joint.diversity.1M <- read.table("maps/rep_1.joint.diversity.block.1.1Mb.0-29999999.bedgraph", header = T)
joint.diversity.1M$avg <- apply(joint.diversity.1M[4:ncol(joint.diversity.1M)], 1, mean)

# loading inferred landscapes
joint.tmrca.1M <- read.table("maps/rep_1.joint.TMRCA.block.1.1Mb.0-29999999.bedgraph", header = T)
joint.tmrca.1M$avg <- apply(joint.tmrca.1M[4:ncol(joint.tmrca.1M)], 1, mean)
joint.rho.1M <- read.table("maps/rep_1.joint.rho.block.1.1Mb.0-29999999.bedgraph", header = T)
joint.theta.1M <- read.table("maps/rep_1.joint.theta.block.1.1Mb.0-29999999.bedgraph", header = T)

inf.lands <- as.data.frame(cbind(joint.diversity.1M$avg, joint.theta.1M$sample_mean, joint.rho.1M$sample_mean, joint.tmrca.1M$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

# does same analyses using the TRUE, SIMULATED landscapes
sim.lands <- as.data.frame(cbind(joint.diversity.1M$avg, sim.theta.1M$sim, sim.rho.1M$sim, sim.tmrca.1M$tmrca))
names(sim.lands) <- c("diversity", "theta", "rho", "tmrca")

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 1000, y = value)) + theme.blank
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, 2 * inf.lands$rho, sim.rho.1M$sim)) # adds simulated mut. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 1000, y = value, colour = variable)) + theme.blank
rho.map <- rho.map + geom_line(data = molten.rho)
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.0025, label = "Cor = 0.92")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.1M$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 1000, y = value, colour = variable)) + theme.blank
theta.map <- theta.map + geom_line(data = molten.theta)
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.006, label = "Cor = 0.94")

tmrca.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$tmrca, sim.tmrca.1M$tmrca))
names(tmrca.plot) <- c("bin", "inf", "sim")
molten.tmrca <- melt(tmrca.plot, id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 1000, y = value, colour = variable)) + theme.blank
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F) + no.legend
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
tmrca.map <- tmrca.map + annotate("text", x = 15000, y = 1.125, label = "Cor = 0.56")


# all together now
theme_set(theme_cowplot(font_size = 12))
legend <- get_legend(rho.map + theme(legend.position="bottom"))
p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
lands.plot.1M <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::save_plot("landscapes.1Mb.phased.pdf", lands.plot.1M, base_width = 9, base_height = 15)




# Linear model 

# centering
inf.lands$thetaC <- inf.lands$theta - mean(inf.lands$theta)
inf.lands$tmrcaC <- inf.lands$tmrca - mean(inf.lands$tmrca)
inf.lands$rhoC <- inf.lands$rho - mean(inf.lands$rho)


plot(diversity~theta, data = inf.lands)
plot(diversity~tmrca, data = inf.lands)
plot(diversity~rho, data = inf.lands)  

plot(theta~tmrca, data = inf.lands)
plot(tmrca~rho, data = inf.lands)
plot(theta~rho, data = inf.lands)

m.diversity <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands)

plot(y=resid(m.diversity, type = "pearson"), x=fitted(m.diversity))
dwtest(m.diversity)
hmctest(m.diversity)
hist(resid(m.diversity))

summary(m.diversity)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.0035935  0.0004619  -7.780 2.97e-08 ***
# theta        0.9616656  0.0385908  24.920  < 2e-16 ***
# rho         -0.1731741  0.1336117  -1.296    0.206    
# tmrca        0.0039649  0.0005119   7.746 3.23e-08 ***
# Multiple R-squared: 0.969,	Adjusted R-squared:  0.9654 

# type 2 anova
anova.diversity <- Anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# thetaC        2.7676e-05  1 1094.3306 0.00000 0.87809
# rhoC          9.7000e-09  1    0.3847 0.54073 0.00031
# tmrcaC        2.6738e-06  1  105.7238 0.00000 0.08483
# thetaC:tmrcaC 5.2650e-07  1   20.8185 0.00012 0.01670
# Residuals     6.3230e-07 25                   0.02006

cor.tab[3,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                  anova.diversity$VarExp[1] * 100, 0, anova.diversity$VarExp[3] * 100, 1000)

g.diversity <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC,
                   data = inf.lands, corr = corAR1(0, ~bin))

summary(g.diversity)
# (Intercept)   0.0033747 0.0000422 79.91382  0.0000
# thetaC        0.9396730 0.0281665 33.36143  0.0000
# rhoC          0.1067068 0.1018576  1.04761  0.3048
# tmrcaC        0.0047202 0.0004235 11.14652  0.0000
# thetaC:tmrcaC 2.1533778 0.4428766  4.86225  0.0001

vif(g.diversity)
# thetaC          rhoC        tmrcaC thetaC:tmrcaC 
# 1.120283      1.303918      1.126846      1.368307

###################################################
#
# Linear model using simulated (true) landscapes
#
###################################################


# 50kb
sim.lands.50kb <- as.data.frame(cbind(joint.diversity.50k$avg, sim.theta.50k$sim, sim.rho.50k$sim, sim.tmrca.50k$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")


hist(sim.lands.50kb$diversity)
hist(sim.lands.50kb$rho)
hist(sim.lands.50kb$theta)
hist(sim.lands.50kb$tmrca)

cor.test(~theta+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.50kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.50kb, method = "spearman") 

plot(theta~tmrca, data = sim.lands.50kb)
plot(theta~rho, data = sim.lands.50kb)
plot(tmrca~rho, data = sim.lands.50kb)

plot(diversity~theta, data = sim.lands.50kb)
plot(diversity~tmrca, data = sim.lands.50kb)
plot(diversity~rho, data = sim.lands.50kb)

# centering
sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)

sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)

m.div.50kb <- lm(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC, data = sim.lands.50kb)
plot(resid(m.div.50kb)~fitted(m.div.50kb))
dwtest(m.div.50kb)
hmctest(m.div.50kb)
hist(resid(m.div.50kb))

summary(m.div.50kb)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   3.346e-03  1.557e-05 214.941   <2e-16 ***
# thetaC        9.656e-01  8.102e-03 119.181   <2e-16 ***
# tmrcaC        3.329e-03  5.719e-05  58.204   <2e-16 ***
# rhoC          3.861e-03  1.005e-02   0.384    0.701    
# thetaC:tmrcaC 1.052e+00  3.489e-02  30.161   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.013678      1.013298      1.006753      1.006017

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

cor.tab[4,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[4]) * 100,
                  anova.diversity$VarExp[1] * 100, 0, anova.diversity$VarExp[2] * 100, 50)

# 
g.diversity.50kb <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                         data = sim.lands.50kb, cor = corAR1(0, ~bin), method = "ML")

summary(g.diversity.50kb)
# Coefficients:
# Value  Std.Error   t-value p-value
# (Intercept)   0.0033459 0.00002474 135.24589  0.0000
# thetaC        0.9497193 0.01092082  86.96412  0.0000
# tmrcaC        0.0033591 0.00005614  59.83852  0.0000
# rhoC          0.0023547 0.01154079   0.20403  0.8384
# thetaC:tmrcaC 1.0541617 0.03239511  32.54076  0.0000


# 200kb
sim.lands.200kb <- as.data.frame(cbind(joint.diversity.200k$avg, sim.theta.200k$sim, sim.rho.200k$sim, sim.tmrca.200k$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

sim.lands.200kb$theta <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
sim.lands.200kb$tmrca <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
sim.lands.200kb$rho <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)

hist(sim.lands.200kb$diversity)
hist(sim.lands.200kb$rho)
hist(sim.lands.200kb$theta)
hist(sim.lands.200kb$tmrca)

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

plot(theta~tmrca, data = sim.lands.200kb)
plot(theta~rho, data = sim.lands.200kb)
plot(tmrca~rho, data = sim.lands.200kb)

plot(diversity~theta, data = sim.lands.200kb)
plot(diversity~tmrca, data = sim.lands.200kb)
plot(diversity~rho, data = sim.lands.200kb)

# centering
sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)

sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)

m.div.200kb <- lm(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   3.364e-03  2.258e-05 148.947   <2e-16 ***
# thetaC        9.774e-01  1.329e-02  73.564   <2e-16 ***
# tmrcaC        3.256e-03  1.232e-04  26.434   <2e-16 ***
# rhoC          6.700e-03  1.949e-02   0.344    0.732    
# thetaC:tmrcaC 1.070e+00  8.639e-02  12.387   <2e-16 ***

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

cor.tab[5,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[4]) * 100,
                  anova.diversity$VarExp[1] * 100, 0, anova.diversity$VarExp[2] * 100, 200)

# 
g.diversity.200kb <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                         data = sim.lands.200kb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.diversity.200kb)
# Coefficients:
# Value  Std.Error   t-value p-value
# (Intercept)   0.0033430 0.00001556 214.82327  0.0000
# thetaC        0.9788779 0.00918391 106.58617  0.0000
# tmrcaC        0.0032101 0.00016997  18.88659  0.0000
# rhoC          0.0116600 0.01312824   0.88816  0.3759
# thetaC:tmrcaC 1.0729407 0.11090158   9.67471  0.0000

vif(g.diversity.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.002986      1.011348      1.016106      1.004828 



# 1Mb
sim.lands.1Mb <- as.data.frame(cbind(joint.diversity.1M$avg, sim.theta.1M$sim, sim.rho.1M$sim, sim.tmrca.1M$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

sim.lands.1Mb$theta <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
sim.lands.1Mb$tmrca <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
sim.lands.1Mb$rho <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)

hist(sim.lands.1Mb$diversity)
hist(sim.lands.1Mb$rho)
hist(sim.lands.1Mb$theta)
hist(sim.lands.1Mb$tmrca)

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

plot(theta~tmrca, data = sim.lands.1Mb)
plot(theta~rho, data = sim.lands.1Mb)
plot(tmrca~rho, data = sim.lands.1Mb)

plot(diversity~theta, data = sim.lands.1Mb)
plot(diversity~tmrca, data = sim.lands.1Mb)
plot(diversity~rho, data = sim.lands.1Mb)

# centering
sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)

sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)

m.div.1Mb <- lm(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   3.402e-03  3.967e-05  85.759  < 2e-16 ***
# thetaC        9.775e-01  3.938e-02  24.821  < 2e-16 ***
# tmrcaC        3.006e-03  4.736e-04   6.347 1.21e-06 ***
# rhoC          5.781e-03  6.109e-02   0.095    0.925    
# thetaC:tmrcaC 2.755e-02  5.086e-01   0.054    0.957    

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

cor.tab[6,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[4]) * 100,
                  anova.diversity$VarExp[1] * 100, 0, anova.diversity$VarExp[2] * 100, 1000)

# 
g.diversity.1Mb <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                       data = sim.lands.1Mb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.diversity.1Mb)
# Coefficients:
# Value  Std.Error   t-value p-value
# (Intercept)    0.0034117 0.00001947 175.20689  0.0000
# thetaC         1.0047440 0.02517057  39.91741  0.0000
# tmrcaC         0.0030554 0.00022456  13.60658  0.0000
# rhoC          -0.0549662 0.04303333  -1.27729  0.2132
# thetaC:tmrcaC -0.0894952 0.24683676  -0.36257  0.7200

vif(g.diversity.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.306371      1.191145      1.202768      1.189138

########################################
#
# R2 Plot
#
########################################

cor.tab.2 <- as.data.frame(cbind(apply(cor.tab, 2, as.numeric)))
names(cor.tab.2)[5] <- "bin.size"
cor.tab.2$type <- c(rep("inf", 3), rep("sim", 3))

molten.cor <- melt(cor.tab.2, id.vars = c("bin.size", "type"))
cor.plot <- ggplot(data = molten.cor, aes(x = bin.size, y = value, colour = variable, fill = type))
cor.plot <- cor.plot + geom_line(data = molten.cor)
cor.plot <- cor.plot + geom_point(aes(shape = type, colour = variable), size = 7)
cor.plot <- cor.plot + scale_x_continuous(breaks = c(50, 200, 1000))
cor.plot <- cor.plot + scale_y_continuous(breaks = pretty_breaks())
cor.plot <- cor.plot + labs(title = NULL, x = "Bin Size (kb)", y = "Variance Explained (%)") + no.legend
cor.plot <- cor.plot + theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))

ggsave("lm.cor.sim.pdf", cor.plot, width = 7, height = 7)

