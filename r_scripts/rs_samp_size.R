# Created: 17/06/2019
# Last modified: 12/08/2020
# Author: Gustavo Barroso


library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)
library(ppls)
library(scales)
library(cowplot)
library(MASS)
library(lmtest)
library(relaimpo)
library(ape)
library(nlme)
library(plyr)

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
# diversity ~ rho + theta + tmrca --- 50 kb scale (UN-PHASED)
#
###################################################

no.x.axis <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
no.legend <- theme(legend.position = "none", legend.title = NULL)

# sim landscapes 50 kb
sim.rho.50k <- read.table("raw_data/rig.sims.rho.50000.bins.txt")
sim.theta.50k <- read.table("raw_data/rig.sims.theta.50000.bins.txt")

# loading computed diversity landscape
unphased.diversity.50k <- read.table("maps/rep_1.unphased.diversity.block.1.50kb.0-29999999.bedgraph", header = T)
unphased.diversity.50k$avg <- apply(unphased.diversity.50k[4:ncol(unphased.diversity.50k)], 1, mean)

# loading inferred landscapes
unphased.tmrca.50k <- read.table("maps/rep_1.unphased.TMRCA.block.1.50kb.0-29999999.bedgraph", header = T)
unphased.tmrca.50k$avg <- apply(unphased.tmrca.50k[4:ncol(unphased.tmrca.50k)], 1, mean)
unphased.rho.50k <- read.table("maps/rep_1.unphased.rho.block.1.50kb.0-29999999.bedgraph", header = T)
unphased.theta.50k <- read.table("maps/rep_1.unphased.theta.block.1.50kb.0-29999999.bedgraph", header = T)

# simulated and inferred maps are highly correlated
cor.test(x = unphased.rho.50k$sample_mean, y = sim.rho.50k$sim, method = "spearman")$estimate
cor.test(x = unphased.theta.50k$sample_mean, y = sim.theta.50k$sim, method = "spearman")$estimate

inf.lands <- as.data.frame(cbind(unphased.diversity.50k$avg, unphased.theta.50k$sample_mean, unphased.rho.50k$sample_mean, unphased.tmrca.50k$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

cor.mat <- cor(inf.lands[,1:4], method = "spearman")
# non-significant correlations
cor.mat[3, 1] <- 0
cor.mat[3, 2] <- 0
cor.mat[4, 3] <- 0
cor.mat[3, 4] <- 0

pdf("sim.cor.50kb.pdf")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 50, y = value, colour = "#F8766D")) 
diversity.map <- diversity.map + geom_line(data = molten.diversity)
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + theme.blank
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi)) + no.legend
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) 

rho.plot <- as.data.frame(cbind(inf.lands$bin, 2 * inf.lands$rho, sim.rho.50k$sim)) # adds simulated rec. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 50, y = value, colour = variable)) 
rho.map <- rho.map + geom_line(data = molten.rho) + scale_fill_manual(values = c("#fc8d59", "#99d594")) + theme.blank
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.009, label = "Cor = 0.92")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.50k$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 50, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.01, label = "Cor = 0.90")

tmrca.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$tmrca, sim.tmrca.50k$tmrca))
names(tmrca.plot) <- c("bin", "inf", "sim")
molten.tmrca <- melt(tmrca.plot, id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 50, y = value, colour = variable)) + theme.blank
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F) + no.legend
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

# all together now
theme_set(theme_cowplot(font_size = 12))
legend <- get_legend(rho.map + theme(legend.position="bottom"))
p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
lands.plot.50k <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::save_plot("landscapes.50kb.unphased.pdf", lands.plot, base_width = 9, base_height = 15)



# Linear model 
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) #***
hist(resid(m.diversity.bc), nclass = 30) # not too bad?
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # ***
Box.test(resid(m.diversity)[order(predict(m.diversity))], type = "Ljung-Box") #***

summary(m.diversity.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.002908   0.001066 -1878.693  < 2e-16 ***
# theta       11.968188   0.120998    98.912  < 2e-16 ***
# rho          1.432329   0.379820     3.771 0.000179 ***
# tmrca        0.088523   0.001348    65.656  < 2e-16 ***
# Multiple R-squared:  0.9701,	Adjusted R-squared:   0.97

anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta       1 0.66540 0.66540 15039.111 0.000000 0.75362
# rho         1 0.00044 0.00044    10.053 0.001599 0.00050
# tmrca       1 0.19073 0.19073  4310.741 0.000000 0.21601
# Residuals 596 0.02637 0.00004                    0.02987



# Because of auto-correlation we compute p-values for the variables using a GLS
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = inf.lands, corr = corAR1(0, ~bin))

summary(g.diversity)
# Parameter estimate(s):
#  Phi 
# 0.5035707 

# Coefficients:
#  Value  Std.Error   t-value p-value
# (Intercept) -0.0041156 0.00011986 -34.33830  0.0000
# theta        0.7410263 0.01755293  42.21667  0.0000
# rho         -0.0710679 0.04499289  -1.57954  0.1147
# tmrca        0.0056709 0.00012527  45.26923  0.0000



###################################################
#
# diversity ~ rho + theta + tmrca --- 200 kb scale (UN-PHASED)
#
###################################################

# sim landscapes 200 kb
sim.rho.200k <- read.table("raw_data/rig.sims.rho.2e+05.bins.txt")
sim.theta.200k <- read.table("raw_data/rig.sims.theta.2e+05.bins.txt")

# loading computed diversity landscape
unphased.diversity.200k <- read.table("maps/rep_1.unphased.diversity.block.1.200kb.0-29999999.bedgraph", header = T)
unphased.diversity.200k$avg <- apply(unphased.diversity.200k[4:ncol(unphased.diversity.200k)], 1, mean)

# loading inferred landscapes
unphased.tmrca.200k <- read.table("maps/rep_1.unphased.TMRCA.block.1.200kb.0-29999999.bedgraph", header = T)
unphased.tmrca.200k$avg <- apply(unphased.tmrca.200k[4:ncol(unphased.tmrca.200k)], 1, mean)
unphased.rho.200k <- read.table("maps/rep_1.unphased.rho.block.1.200kb.0-29999999.bedgraph", header = T)
unphased.theta.200k <- read.table("maps/rep_1.unphased.theta.block.1.200kb.0-29999999.bedgraph", header = T)

# simulated and inferred maps are highly correlated
cor.test(x = unphased.rho.200k$sample_mean, y = sim.rho.200k$sim, method = "spearman")$estimate
cor.test(x = unphased.theta.200k$sample_mean, y = sim.theta.200k$sim, method = "spearman")$estimate

inf.lands <- as.data.frame(cbind(unphased.diversity.200k$avg, unphased.theta.200k$sample_mean, unphased.rho.200k$sample_mean, unphased.tmrca.200k$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

cor.mat <- cor(inf.lands[,1:4], method = "spearman")$estimate
# non-significant correlations
cor.mat[3, 1] <- 0
cor.mat[3, 2] <- 0
cor.mat[4, 3] <- 0
cor.mat[3, 4] <- 0

pdf("sim.cor.200kb.pdf")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 200, y = value, colour = "#F8766D")) 
diversity.map <- diversity.map + geom_line(data = molten.diversity)
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$rho, sim.rho.200k$sim)) # adds simulated rec. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 200, y = value, colour = variable)) 
rho.map <- rho.map + geom_line(data = molten.rho)
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.0075, label = "Cor = 0.92")

theta.plot <-  as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.200k$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 200, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta)
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.01, label = "Cor = 0.93")

tmrca.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$tmrca, sim.tmrca.200k$tmrca))
names(tmrca.plot) <- c("bin", "inf", "sim")
molten.tmrca <- melt(tmrca.plot, id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 50, y = value, colour = variable)) 
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F) + no.legend
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

# all together now
theme_set(theme_cowplot(font_size = 12))
legend <- get_legend(rho.map + theme(legend.position="bottom"))
p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
lands.plot.200k <-  plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::ggsave("landscapes.200kb.unphased.pdf", lands.plot, width = 9, height = 15)




# Linear model 

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands)
plot(m.diversity, which = 2)

shapiro.test(resid(m.diversity)) #**
hist(resid(m.diversity), nclass = 30) # not too bad?
hmctest(m.diversity, nsim = 3000) # NS
dwtest(m.diversity) # **
Box.test(resid(m.diversity)[order(predict(m.diversity))], type = "Ljung-Box") #**

summary(m.diversity)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.0039790  0.0001911 -20.824   <2e-16 ***
# theta        0.7315314  0.0180289  40.575   <2e-16 ***
# rho         -0.0564721  0.0622994  -0.906    0.366    
# tmrca        0.0055288  0.0002552  21.662   <2e-16 ***
# Multiple R-squared:  0.9514,	Adjusted R-squared:  0.9504

anova.diversity <- Anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# theta     0.00033926   1 1646.3666 0.00000 0.72770
# rho       0.00000017   1    0.8217 0.36618 0.00036
# tmrca     0.00009669   1  469.2271 0.00000 0.20740
# Residuals 0.00003009 146                   0.06453


# Because of auto-correlation we compute p-values for the variables using a GLS
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = inf.lands, corr = corAR1(0, ~bin))

summary(g.diversity)
# Parameter estimate(s):
#  Phi 
# 0.2528086 

# Coefficients:
#  Value  Std.Error   t-value p-value
# (Intercept) -0.0040086 0.00019120 -20.96514  0.0000
# theta        0.7339877 0.02100331  34.94629  0.0000
# rho         -0.0917765 0.06397296  -1.43461  0.1535
# tmrca        0.0055893 0.00025527  21.89586  0.0000


###################################################
#
# diversity ~ rho + theta + tmrca --- 1Mb scale (UN-PHASED)
#
###################################################

# sim landscapes
sim.rho.1M <- read.table("raw_data/rig.sims.rho.1e+06.bins.txt")
sim.theta.1M <- read.table("raw_data/rig.sims.theta.1e+06.bins.txt")

# loading computed diversity landscape
unphased.diversity.1M <- read.table("maps/rep_1.unphased.diversity.block.1.1Mb.0-29999999.bedgraph", header = T)
unphased.diversity.1M$avg <- apply(unphased.diversity.1M[4:ncol(unphased.diversity.1M)], 1, mean)

# loading inferred landscapes
unphased.tmrca.1M <- read.table("maps/rep_1.unphased.TMRCA.block.1.1Mb.0-29999999.bedgraph", header = T)
unphased.tmrca.1M$avg <- apply(unphased.tmrca.1M[4:ncol(unphased.tmrca.1M)], 1, mean)
unphased.rho.1M <- read.table("maps/rep_1.unphased.rho.block.1.1Mb.0-29999999.bedgraph", header = T)
unphased.theta.1M <- read.table("maps/rep_1.unphased.theta.block.1.1Mb.0-29999999.bedgraph", header = T)

# simulated and inferred maps are highly correlated
cor.test(x = unphased.rho.1M$sample_mean, y = sim.rho.1M$sim, method = "spearman")$estimate
cor.test(x = unphased.theta.1M$sample_mean, y = sim.theta.1M$sim, method = "spearman")$estimate

inf.lands <- as.data.frame(cbind(unphased.diversity.1M$avg, unphased.theta.1M$sample_mean, unphased.rho.1M$sample_mean, unphased.tmrca.1M$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

cor.mat <- cor(inf.lands[,1:4], method = "spearman")$estimate
# non-significant correlations
cor.mat[3, 1] <- 0
cor.mat[3, 2] <- 0
cor.mat[4, 3] <- 0
cor.mat[3, 4] <- 0

pdf("sim.cor.1Mb.pdf")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 1000, y = value, colour = "#F8766D")) 
diversity.map <- diversity.map + geom_line(data = molten.diversity)
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$rho, sim.rho.1M$sim)) # adds simulated mut. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 1000, y = value, colour = variable)) 
rho.map <- rho.map + geom_line(data = molten.rho)
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.0025, label = "Cor = 0.90")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.1M$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 1000, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta)
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.0085, label = "Cor = 0.93")

tmrca.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$tmrca, sim.tmrca.1M$tmrca))
names(tmrca.plot) <- c("bin", "inf", "sim")
molten.tmrca <- melt(tmrca.plot, id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 50, y = value, colour = variable)) 
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F) + no.legend
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

# all together now
theme_set(theme_cowplot(font_size = 12))
legend <- get_legend(rho.map + theme(legend.position="bottom"))
p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
lands.plot.1M <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::ggsave("landscapes.1Mb.unphased.pdf", lands.plot, width = 9, height = 15)






# Linear model WITHOUT interactions

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(-0.0, 1., len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.int.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # NS
hist(resid(m.diversity.bc), nclass = 30) 
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # NS
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.006445   0.009025 -222.319  < 2e-16 ***
# theta       11.846176   0.529758   22.361  < 2e-16 ***
# rho          0.642341   2.159814    0.297    0.769    
# tmrca        0.086390   0.011601    7.447  6.6e-08 ***
# Multiple R-squared:   0.9597,	Adjusted R-squared:  0.9551

anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df     Sum Sq    Mean Sq  F value    Pr(>F)  VarExp
# theta      1 0.0097537 0.0097537 557.6809 0.0000000 0.86412
# rho        1 0.0001091 0.0001091   6.2385 0.0191573 0.00967
# tmrca      1 0.0009699 0.0009699  55.4557 0.0000001 0.08593
# Residuals 26 0.0004547 0.0000175                    0.04029

relimp.diversity <- calc.relimp(m.diversity.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 0.0003892234 
# Proportion of variance explained by model: 95.97%
# Relative importance metrics: 
# theta 0.86470529
# rho   0.02175124
# tmrca 0.11354347






# Linear model with interactions

m.diversity.int <- lm(diversity ~ (theta + rho + tmrca) ^ 2, data = inf.lands)
plot(m.diversity.int, which = 2)

bc.diversity <- boxcox(m.diversity.int, lambda = seq(0.5, 1.5, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.int.bc <- update(m.diversity.int, (diversity^l -1)/l~.)
plot(m.diversity.int.bc, which = 2)

shapiro.test(resid(m.diversity.int.bc)) #*
hist(resid(m.diversity.int.bc), nclass = 30)
hmctest(m.diversity.int.bc, nsim = 3000) # NS
dwtest(m.diversity.int.bc) # NS
Box.test(resid(m.diversity.int.bc)[order(predict(m.diversity.int.bc))], type = "Ljung-Box") # NS

summary(m.diversity.int.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -1.129078   0.003216 -351.066  < 2e-16 ***
# theta        -0.567943   0.598751   -0.949  0.35271    
# rho          -1.560005   3.876741   -0.402  0.69110    
# tmrca        -0.004208   0.003848   -1.094  0.28546    
# theta:rho    23.703281 200.264165    0.118  0.90681    
# theta:tmrca   2.702974   0.752503    3.592  0.00154 ** 
# rho:tmrca     2.531711   4.458653    0.568  0.57566    
# Multiple R-squared:  0.9792,	Adjusted R-squared:  0.9738 

anova.diversity <- anova(m.diversity.int.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta        1 1.3125e-04 1.3125e-04 959.9745 0.00000 0.86700
# rho          1 2.0020e-06 2.0020e-06  14.6455 0.00086 0.01323
# tmrca        1 1.1880e-05 1.1880e-05  86.8881 0.00000 0.07847
# theta:rho    1 1.1420e-06 1.1420e-06   8.3542 0.00825 0.00755
# theta:tmrca  1 1.9210e-06 1.9210e-06  14.0464 0.00105 0.01269
# rho:tmrca    1 4.4000e-08 4.4000e-08   0.3224 0.57566 0.00029
# Residuals   23 3.1450e-06 1.3700e-07                  0.02077

relimp.diversity <- calc.relimp(m.diversity.int.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 5.220301e-06 
# Proportion of variance explained by model: 97.92%
# Relative importance metrics: 
# theta       0.855793450
# rho         0.018899581
# tmrca       0.104145431
# theta:rho   0.003578760
# theta:tmrca 0.014876883
# rho:tmrca   0.002705896


# Comparing the model with interactions and the model without
anova(m.diversity.bc, m.diversity.int.bc) #*


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

cor.mat <- cor(inf.lands[,1:4], method = "spearman")
# non-significant correlations
cor.mat[3, 1] <- 0
cor.mat[3, 2] <- 0
cor.mat[4, 3] <- 0


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

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) #***
hist(resid(m.diversity.bc), nclass = 30) # not too bad?
hmctest(m.diversity.bc, nsim = 3000) # *
dwtest(m.diversity.bc) # ***
Box.test(resid(m.diversity)[order(predict(m.diversity))], type = "Ljung-Box") #***


summary(m.diversity.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.9309400  0.0008377 -2305.14  < 2e-16 ***
# theta       13.9738081  0.1280305   109.14  < 2e-16 ***
# rho          1.1733146  0.3760272     3.12  0.00189 ** 
# tmrca        0.0653199  0.0008903    73.37  < 2e-16 ***
# Multiple R-squared:  0.9762,	Adjusted R-squared:  0.9761

# type 2 anova
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)
anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta     0.300157   1 11912.4867 0.0000000 0.66546
# rho       0.000245   1     9.7362 0.0018941 0.00054
# tmrca     0.135632   1  5382.8809 0.0000000 0.30070
# Residuals 0.015017 596                      0.03329

cor.tab[1,] <- c(66.5 + 30, 66.5, 0, 30, 50)

# Because of auto-correlation we compute p-values for the variables using a GLS
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = inf.lands, corr = corAR1(0, ~bin))

summary(g.diversity)
# Parameter estimate(s):
#  Phi 
# 0.2579507 

# Coefficients:
#  Value  Std.Error   t-value p-value
# (Intercept) -0.0037347 0.00049822 -7.496035  0.0000
# theta        0.9737686 0.03737108 26.056741  0.0000
# rho         -0.1307381 0.12647873 -1.033676  0.3108
# tmrca        0.0040447 0.00053973  7.493863  0.0000



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

cor.mat <- cor(inf.lands[,1:4], method = "spearman")
# non-significant correlations
cor.mat[3, 1] <- 0
cor.mat[3, 2] <- 0
cor.mat[4, 3] <- 0
cor.mat[3, 4] <- 0

# does same analyses using the TRUE, SIMULATED landscapes
sim.lands <- as.data.frame(cbind(joint.diversity.200k$avg, sim.theta.200k$sim, sim.rho.200k$sim, sim.tmrca.200k$tmrca))
names(sim.lands) <- c("diversity", "theta", "rho", "tmrca")
cor.mat.truth <- cor(sim.lands[,1:4], method = "spearman")

# adds it to upper triangle of cor.mat
cor.mat[1, 2] <- 0.93
cor.mat[1, 3] <- 0
cor.mat[1, 4] <- 0.36
cor.mat[2, 3] <- 0
cor.mat[2, 4] <- 0
cor.mat[3, 4] <- 0

pdf("sim.cor.200kb.pdf")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "full")
dev.off()

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
cowplot::save_plot("landscapes.200kb.phased.pdf", lands.plot, base_width = 9, base_height = 15)



# Linear model 

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0., 1., len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.int.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) #**
hist(resid(m.diversity.bc), nclass = 30) # not too bad?
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # **
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.785113   0.001346 -1326.384   <2e-16 ***
# theta       11.181677   0.171113    65.347   <2e-16 ***
# rho         -0.088707   0.559844    -0.158    0.874    
# tmrca        0.048154   0.001500    32.102   <2e-16 *
# Multiple R-squared:  0.9801,	Adjusted R-squared:  0.9797

# type 2 anova
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)
anova.diversity
# theta     0.041214   1 4270.2126 0.00000 0.78399
# rho       0.000000   1    0.0251 0.87432 0.00000
# tmrca     0.009946   1 1030.5165 0.00000 0.18920
# Residuals 0.001409 146                   0.02680

cor.tab[2,] <- c(78.4 + 18.9, 78.4, 0, 18.9, 200)

# Because of auto-correlation we compute p-values for the variables using a GLS
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = inf.lands, corr = corAR1(0, ~bin))

summary(g.diversity)
# Parameter estimate(s):
#  Phi 
# 0.2078107 

# Coefficients:
#  Value  Std.Error   t-value p-value
# (Intercept) -0.0037655 0.00016568 -22.72800  0.0000
# theta        0.9731124 0.02373102  41.00593  0.0000
# rho         -0.1274976 0.07091808  -1.79782  0.0743
# tmrca        0.0040808 0.00018333  22.25969  0.0000


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

cor.mat <- cor(inf.lands[,1:4], method = "spearman")
# non-significant correlations
cor.mat[3, 1] <- 0
cor.mat[3, 2] <- 0
cor.mat[4, 3] <- 0
cor.mat[3, 4] <- 0

# does same analyses using the TRUE, SIMULATED landscapes
sim.lands <- as.data.frame(cbind(joint.diversity.1M$avg, sim.theta.1M$sim, sim.rho.1M$sim, sim.tmrca.1M$tmrca))
names(sim.lands) <- c("diversity", "theta", "rho", "tmrca")
cor.mat.truth <- cor(sim.lands[,1:4], method = "spearman")

# adds it to upper triangle of cor.mat
cor.mat[1, 2] <- 0.95
cor.mat[1, 3] <- 0
cor.mat[1, 4] <- 0.33
cor.mat[2, 3] <- 0
cor.mat[2, 4] <- 0
cor.mat[3, 4] <- 0

pdf("sim.cor.1Mb.pdf")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "full")
dev.off()

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
cowplot::save_plot("landscapes.1Mb.phased.pdf", lands.plot, base_width = 9, base_height = 15)




# Linear model 

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands)
plot(m.diversity, which = 2)

shapiro.test(resid(m.diversity)) # NS
hist(resid(m.diversity.bc), nclass = 30) 
hmctest(m.diversity, nsim = 3000) # NS
dwtest(m.diversity) # NS
Box.test(resid(m.diversity.bc)[order(predict(m.diversity))], type = "Ljung-Box") # *

summary(m.diversity)
# Coefficients:
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
# theta     2.7676e-05  1 620.9843 0.00000 0.87628
# rho       7.4900e-08  1   1.6799 0.20633 0.00237
# tmrca     2.6738e-06  1  59.9936 0.00000 0.08466
# Residuals 1.1588e-06 26                  0.03669

cor.tab[3,] <- c(87.6 + 8.5, 87.6, 0, 8.5, 1000)




########################################
#
# Lands Plot
#
#######################################

lands.plot.all <- plot_grid(lands.plot.50k, lands.plot.200k, lands.plot.1M, ncol = 3, labels = "AUTO", label_size = 36)
cowplot::save_plot("sim.landscapes.pdf", lands.plot.all, base_width = 27, base_height = 15)

###################################################
#
# Linear model using simulated (true) landscapes
#
###################################################

# 50kb
sim.lands <- as.data.frame(cbind(joint.diversity.50k$avg, sim.theta.50k$sim, sim.rho.50k$sim, sim.tmrca.50k$tmrca))
names(sim.lands) <- c("diversity", "theta", "rho", "tmrca")

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = sim.lands)
bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)
summary(m.diversity.bc)
# Coefficients:
# (Intercept) -1.5418649  0.0006435 -2395.998   <2e-16 ***
# theta        6.7242236  0.0826287    81.379   <2e-16 ***
# rho          0.1053594  0.1023716     1.029    0.304    
# tmrca        0.0233221  0.0005827    40.027   <2e-16 ***
# Multiple R-squared:  0.938,	Adjusted R-squared:  0.9377 
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)
anova.diversity
# theta     0.098780   1 6622.5077 0.00000 0.75070
# rho       0.000016   1    1.0592 0.30381 0.00012
# tmrca     0.023898   1 1602.1907 0.00000 0.18162
# Residuals 0.008890 596                   0.06756

cor.tab[4,] <- c(75.1 + 18.2, 75.1, 0, 18.2, 50)



# 200kb
sim.lands <- as.data.frame(cbind(joint.diversity.200k$avg, sim.theta.200k$sim, sim.rho.200k$sim, sim.tmrca.200k$tmrca))
names(sim.lands) <- c("diversity", "theta", "rho", "tmrca")

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = sim.lands)
bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)
summary(m.diversity.bc)
# Coefficients:
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept) -1.3222174  0.0006662 -1984.62   <2e-16 ***
# theta        3.7357895  0.0701360    53.27   <2e-16 ***
# rho         -0.0294771  0.1017289    -0.29    0.772    
# tmrca        0.0127860  0.0006481    19.73   <2e-16 ***
# Multiple R-squared:  0.9602,	Adjusted R-squared:  0.9594 
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)
anova.diversity
# theta     0.0059432   1 2837.152 0.00000 0.84129
# rho       0.0000002   1    0.084 0.77241 0.00002
# tmrca     0.0008152   1  389.166 0.00000 0.11540
# Residuals 0.0003058 146                  0.04329

cor.tab[5,] <- c(84.1 + 11.5, 84.1, 0, 11.5, 200)





# 1Mb
sim.lands <- as.data.frame(cbind(joint.diversity.1M$avg, sim.theta.1M$sim, sim.rho.1M$sim, sim.tmrca.1M$tmrca))
names(sim.lands) <- c("diversity", "theta", "rho", "tmrca")

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = sim.lands)
bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)
summary(m.diversity.bc)
# Coefficients:
# Estimate Std. Error   t value Pr(>|t|)    
# (Intercept) -1.0390655  0.0005530 -1878.884  < 2e-16 ***
# theta        1.1869697  0.0453478    26.175  < 2e-16 ***
# rho          0.0097638  0.0715858     0.136    0.893    
# tmrca        0.0036510  0.0005649     6.463 7.54e-07 ***
# Multiple R-squared:  0.9695,	Adjusted R-squared:  0.966 
anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)
anova.diversity
# theta      1 5.0617e-05 5.0617e-05 783.321  0.000 0.91932
# rho        1 6.3000e-08 6.3000e-08   0.969  0.334 0.00114
# tmrca      1 2.6990e-06 2.6990e-06  41.772  0.000 0.04902
# Residuals 26 1.6800e-06 6.5000e-08                0.03051

cor.tab[6,] <- c(91.9 + 4.9, 91.9, 0, 4.9, 1000)


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
cor.plot <- cor.plot + labs(title = NULL, x = "Bin Size (kb)", y = "Var. Explained")
cor.plot <- cor.plot + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

ggsave("lm.cor.sim.pdf", cor.plot, width = 7, height = 7)

