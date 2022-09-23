# Created: 25/10/2020
# Last modified: 29/10/2020
# Author: Gustavo Barroso 
# This script produces the R^2 plots from different simulated scenarios

library(tidyverse)
library(car)
library(scales)
library(cowplot)

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

################################################
#
# bottleneck_r_1e-8_flat_mu
#
################################################

d <- "other_neutral_scenarios/bottleneck_r_1e-8_flat_mu/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.50kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "/sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "/rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "/rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)

  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "rho", "tmrca")

  # centering
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)

  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)

  m.div.50kb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.50kb)

  summary(m.div.50kb)

  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)

  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.200kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "/sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "/rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "/rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.1Mb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "/sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


r2.sim.avg_1 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_1) <- c("Total", "Rho", "TMRCA")
r2.sim.avg_1$Theta <- NA
r2.sim.avg_1 <- r2.sim.avg_1[,c(1,4,2,3)]
r2.sim.avg_1$bin.size <- c(50, 200, 1000)
r2.sim.avg_1$rho <- "r = 1e-8"
r2.sim.avg_1$mu_block <- "mu = flat"

################################################
#
# bottleneck_r_1e-8_mu_change_50kb
#
################################################

d <- "other_neutral_scenarios/bottleneck_r_1e-8_mu_change_50kb/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "/sim.rho.50000.map", sep = ""), header = T)
sim.theta.50kb <- read.table(paste(d, "/sim.theta.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "/sim.rho.2e+05.map", sep = ""), header = T)
sim.theta.200kb <- read.table(paste(d, "/sim.theta.2e+05.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "/rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "/rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.1e+06.map", sep = ""), header = T)
sim.theta.1Mb <- read.table(paste(d,"sim.theta.1e+06.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_2 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_2) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.avg_2$bin.size <- c(50, 200, 1000)
r2.sim.avg_2$rho <- "r = 1e-8"
r2.sim.avg_2$mu_block <- "mu block ~50 kb"


################################################
#
# bottleneck_r_1e-8_mu_change_500kb
#
################################################

d <- "other_neutral_scenarios/bottleneck_r_1e-8_mu_change_500kb/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)
sim.theta.50kb <- read.table(paste(d, "sim.theta.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.2e+05.map", sep = ""), header = T)
sim.theta.200kb <- read.table(paste(d,"sim.theta.2e+05.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.1e+06.map", sep = ""), header = T)
sim.theta.1Mb <- read.table(paste(d,"sim.theta.1e+06.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_3 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_3) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.avg_3$bin.size <- c(50, 200, 1000)
r2.sim.avg_3$rho <- "r = 1e-8"
r2.sim.avg_3$mu_block <- "mu block ~500 kb"


################################################
#
# bottleneck_r_1e-9_flat_mu
#
################################################

d <- "other_neutral_scenarios/bottleneck_r_1e-9_flat_mu/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.50kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.200kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.1Mb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_4 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_4) <- c("Total", "Rho", "TMRCA")
r2.sim.avg_4$Theta <- NA
r2.sim.avg_4 <- r2.sim.avg_4[,c(1,4,2,3)]
r2.sim.avg_4$bin.size <- c(50, 200, 1000)
r2.sim.avg_4$rho <- "r = 1e-9"
r2.sim.avg_4$mu_block <- "mu = flat"


################################################
#
# bottleneck_r_1e-9_mu_change_50kb
#
################################################

d <- "other_neutral_scenarios/bottleneck_r_1e-9_mu_change_50kb/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)
sim.theta.50kb <- read.table(paste(d, "sim.theta.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.2e+05.map", sep = ""), header = T)
sim.theta.200kb <- read.table(paste(d,"sim.theta.2e+05.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.1e+06.map", sep = ""), header = T)
sim.theta.1Mb <- read.table(paste(d,"sim.theta.1e+06.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_5 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_5) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.avg_5$bin.size <- c(50, 200, 1000)
r2.sim.avg_5$rho <- "r = 1e-9"
r2.sim.avg_5$mu_block <- "mu block ~50 kb"

################################################
#
# bottleneck_r_1e-9_mu_change_500kb
#
################################################

d <- "other_neutral_scenarios/bottleneck_r_1e-9_mu_change_500kb/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)
sim.theta.50kb <- read.table(paste(d, "sim.theta.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.2e+05.map", sep = ""), header = T)
sim.theta.200kb <- read.table(paste(d,"sim.theta.2e+05.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.1e+06.map", sep = ""), header = T)
sim.theta.1Mb <- read.table(paste(d,"sim.theta.1e+06.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_6 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_6) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.avg_6$bin.size <- c(50, 200, 1000)
r2.sim.avg_6$rho <- "r = 1e-9"
r2.sim.avg_6$mu_block <- "mu block ~500 kb"

################################################
#
# bottleneck plot
#
################################################

r2.sim.bottleneck <- rbind.data.frame(r2.sim.avg_1, r2.sim.avg_2, r2.sim.avg_3, r2.sim.avg_4, r2.sim.avg_5, r2.sim.avg_6, make.row.names = F)
r2.sim.bottleneck$mu_block = factor(r2.sim.bottleneck$mu_block, levels=c('mu = flat','mu block ~500 kb','mu block ~50 kb'))

molten.r2 <- melt(r2.sim.bottleneck, id.vars = c("bin.size", "rho", "mu_block"))

r2.plot_bottleneck <- ggplot(data = molten.r2, aes(x = bin.size, y = value, shape = variable))
r2.plot_bottleneck <- r2.plot_bottleneck + geom_line(data = molten.r2) + facet_grid(rho ~ mu_block) 
r2.plot_bottleneck <- r2.plot_bottleneck + geom_point(aes(shape = variable), size = 4)
r2.plot_bottleneck <- r2.plot_bottleneck + scale_shape_manual(values = c(0, 1, 2, 8))
r2.plot_bottleneck <- r2.plot_bottleneck + scale_x_continuous(breaks = c(50, 200, 1000), limits = c(40, 1200), trans="log10") 
r2.plot_bottleneck <- r2.plot_bottleneck + scale_y_continuous(breaks = pretty_breaks())
r2.plot_bottleneck <- r2.plot_bottleneck + labs(title = NULL, x = "Bin Size (kb)", y = "Variance Explained (%)") + theme_bw()
r2.plot_bottleneck <- r2.plot_bottleneck + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12))
ggsave("other_neutral_scenarios/r2.bottleneck.pdf", r2.plot_bottleneck, device = "pdf", width = 7, height = 7)

################################################
#
# flat_Ne_r_1e-8_flat_mu
#
################################################

d <- "other_neutral_scenarios/flat_Ne_r_1e-8_flat_mu/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.50kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.200kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.1Mb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_1 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_1) <- c("Total", "Rho", "TMRCA")
r2.sim.avg_1$Theta <- NA
r2.sim.avg_1 <- r2.sim.avg_1[,c(1,4,2,3)]
r2.sim.avg_1$bin.size <- c(50, 200, 1000)
r2.sim.avg_1$rho <- "r = 1e-8"
r2.sim.avg_1$mu_block <- "mu = flat"

################################################
#
# flat_Ne_r_1e-8_mu_change_50kb
#
################################################

d <- "other_neutral_scenarios/flat_Ne_r_1e-8_mu_change_50kb/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)
sim.theta.50kb <- read.table(paste(d, "sim.theta.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.2e+05.map", sep = ""), header = T)
sim.theta.200kb <- read.table(paste(d,"sim.theta.2e+05.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.1e+06.map", sep = ""), header = T)
sim.theta.1Mb <- read.table(paste(d,"sim.theta.1e+06.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_2 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_2) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.avg_2$bin.size <- c(50, 200, 1000)
r2.sim.avg_2$rho <- "r = 1e-8"
r2.sim.avg_2$mu_block <- "mu block ~50 kb"


################################################
#
# flat_Ne_r_1e-8_mu_change_500kb
#
################################################

d <- "other_neutral_scenarios/flat_Ne_r_1e-8_mu_change_500kb/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)
sim.theta.50kb <- read.table(paste(d, "sim.theta.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.2e+05.map", sep = ""), header = T)
sim.theta.200kb <- read.table(paste(d,"sim.theta.2e+05.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.1e+06.map", sep = ""), header = T)
sim.theta.1Mb <- read.table(paste(d,"sim.theta.1e+06.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_3 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_3) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.avg_3$bin.size <- c(50, 200, 1000)
r2.sim.avg_3$rho <- "r = 1e-8"
r2.sim.avg_3$mu_block <- "mu block ~500 kb"


################################################
#
# flat_Ne_r_1e-9_flat_mu
#
################################################

d <- "other_neutral_scenarios/flat_Ne_r_1e-9_flat_mu/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.50kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.200kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.1Mb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ rhoC + tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_4 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_4) <- c("Total", "Rho", "TMRCA")
r2.sim.avg_4$Theta <- NA
r2.sim.avg_4 <- r2.sim.avg_4[,c(1,4,2,3)]
r2.sim.avg_4$bin.size <- c(50, 200, 1000)
r2.sim.avg_4$rho <- "r = 1e-9"
r2.sim.avg_4$mu_block <- "mu = flat"

################################################
#
# flat_Ne_r_1e-9_mu_change_50kb
#
################################################

d <- "other_neutral_scenarios/flat_Ne_r_1e-9_mu_change_50kb/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)
sim.theta.50kb <- read.table(paste(d, "sim.theta.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.2e+05.map", sep = ""), header = T)
sim.theta.200kb <- read.table(paste(d,"sim.theta.2e+05.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.1e+06.map", sep = ""), header = T)
sim.theta.1Mb <- read.table(paste(d,"sim.theta.1e+06.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_5 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_5) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.avg_5$bin.size <- c(50, 200, 1000)
r2.sim.avg_5$rho <- "r = 1e-9"
r2.sim.avg_5$mu_block <- "mu block ~50 kb"

################################################
#
# flat_Ne_r_1e-9_mu_change_500kb
#
################################################

d <- "other_neutral_scenarios/flat_Ne_r_1e-9_mu_change_500kb/"

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table(paste(d, "sim.rho.50000.map", sep = ""), header = T)
sim.theta.50kb <- read.table(paste(d, "sim.theta.50000.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
  sim.lands.50kb <- as.data.frame(cbind(rep_pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_tmrca.50kb$tmrca))
  names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.50kb$thetaC <- sim.lands.50kb$theta - mean(sim.lands.50kb$theta)
  sim.lands.50kb$tmrcaC <- sim.lands.50kb$tmrca - mean(sim.lands.50kb$tmrca)
  sim.lands.50kb$rhoC <- sim.lands.50kb$rho - mean(sim.lands.50kb$rho)
  
  sim.lands.50kb$bin <- 1:nrow(sim.lands.50kb)
  
  m.div.50kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
  
  summary(m.div.50kb)
  
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.50kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.50kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.50kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.50kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))


# 200kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table(paste(d, "sim.rho.2e+05.map", sep = ""), header = T)
sim.theta.200kb <- read.table(paste(d,"sim.theta.2e+05.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste(d, "rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
  sim.lands.200kb <- as.data.frame(cbind(rep_pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_tmrca.200kb$tmrca))
  names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.200kb$thetaC <- sim.lands.200kb$theta - mean(sim.lands.200kb$theta)
  sim.lands.200kb$tmrcaC <- sim.lands.200kb$tmrca - mean(sim.lands.200kb$tmrca)
  sim.lands.200kb$rhoC <- sim.lands.200kb$rho - mean(sim.lands.200kb$rho)
  
  sim.lands.200kb$bin <- 1:nrow(sim.lands.200kb)
  
  m.div.200kb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
  
  summary(m.div.200kb)
  
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.200kb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.200kb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.200kb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.200kb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))


# 1Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table(paste(d, "sim.rho.1e+06.map", sep = ""), header = T)
sim.theta.1Mb <- read.table(paste(d,"sim.theta.1e+06.map", sep = ""), header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste(d, "rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste(d, "rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
  sim.lands.1Mb <- as.data.frame(cbind(rep_pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_tmrca.1Mb$tmrca))
  names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
  
  # centering
  sim.lands.1Mb$thetaC <- sim.lands.1Mb$theta - mean(sim.lands.1Mb$theta)
  sim.lands.1Mb$tmrcaC <- sim.lands.1Mb$tmrca - mean(sim.lands.1Mb$tmrca)
  sim.lands.1Mb$rhoC <- sim.lands.1Mb$rho - mean(sim.lands.1Mb$rho)
  
  sim.lands.1Mb$bin <- 1:nrow(sim.lands.1Mb)
  
  m.div.1Mb <- lm( I(diversity * 1e+6) ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
  
  summary(m.div.1Mb)
  
  anova.diversity <- Anova(m.div.1Mb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.sim.1Mb[1, i] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
  r2.sim.1Mb[2, i] <- anova.diversity$VarExp[1] * 100
  r2.sim.1Mb[3, i] <- anova.diversity$VarExp[2] * 100
  r2.sim.1Mb[4, i] <- anova.diversity$VarExp[3] * 100
}

r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))


# true landscapes
r2.sim.avg_6 <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg_6) <- c("Total", "Theta", "Rho", "TMRCA")
r2.sim.avg_6$bin.size <- c(50, 200, 1000)
r2.sim.avg_6$rho <- "r = 1e-9"
r2.sim.avg_6$mu_block <- "mu block ~500 kb"

################################################
#
# flat_Ne plot
#
################################################

r2.sim.flat_Ne <- rbind.data.frame(r2.sim.avg_1, r2.sim.avg_2, r2.sim.avg_3, r2.sim.avg_4, r2.sim.avg_5, r2.sim.avg_6, make.row.names = F)
r2.sim.flat_Ne$mu_block = factor(r2.sim.flat_Ne$mu_block, levels=c('mu = flat','mu block ~500 kb','mu block ~50 kb'))

molten.r2 <- melt(r2.sim.flat_Ne, id.vars = c("bin.size", "rho", "mu_block"))

r2.plot_flat_Ne <- ggplot(data = molten.r2, aes(x = bin.size, y = value, shape = variable))
r2.plot_flat_Ne <- r2.plot_flat_Ne + geom_line(data = molten.r2) + facet_grid(rho ~ mu_block) 
r2.plot_flat_Ne <- r2.plot_flat_Ne + geom_point(aes(shape = variable), size = 4)
r2.plot_flat_Ne <- r2.plot_flat_Ne + scale_shape_manual(values = c(0, 1, 2, 8))
r2.plot_flat_Ne <- r2.plot_flat_Ne + scale_x_continuous(breaks = c(50, 200, 1000), limits = c(40, 1200), trans="log10") 
r2.plot_flat_Ne <- r2.plot_flat_Ne + scale_y_continuous(breaks = pretty_breaks())
r2.plot_flat_Ne <- r2.plot_flat_Ne + labs(title = NULL, x = "Bin Size (kb)", y = "Variance Explained (%)") + theme_bw()
r2.plot_flat_Ne <- r2.plot_flat_Ne + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12))
ggsave("other_neutral_scenarios/r2.flat_Ne.pdf", r2.plot_flat_Ne, device = "pdf", width = 7, height = 7)


################################################
#
# SLiM simulations with BGS
#
################################################

reps <- 1:10
nreps <- length(reps)

#Regular Region

# true landscapes, mu/r = 10

## 50 kb scale

r2.regular.50kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.regular.50kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.regular.50kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/RegularRegion/1r/RealLandscapes/RealLandscapes_50kb_rep", i, ".csv", sep = "")
  
  sim.50kb <- read.table(path, sep = ",", header = T)
  
  sim.50kb$thetaC <- sim.50kb$MutRate - mean(sim.50kb$MutRate)
  sim.50kb$tmrcaC <- sim.50kb$AverageTmrca - mean(sim.50kb$AverageTmrca)
  m.div.50kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.50kb)
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.regular.50kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.regular.50kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.regular.50kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.regular.50kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

## 200 kb scale

r2.regular.200kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.regular.200kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.regular.200kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/RegularRegion/1r/RealLandscapes/RealLandscapes_200kb_rep", i, ".csv", sep = "")
  
  sim.200kb <- read.table(path, sep = ",", header = T)
  
  sim.200kb$thetaC <- sim.200kb$MutRate - mean(sim.200kb$MutRate)
  sim.200kb$tmrcaC <- sim.200kb$AverageTmrca - mean(sim.200kb$AverageTmrca)
  m.div.200kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.200kb)
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.regular.200kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.regular.200kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.regular.200kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.regular.200kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

## 1000 kb scale

r2.regular.1000kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.regular.1000kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.regular.1000kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/RegularRegion/1r/RealLandscapes/RealLandscapes_1000kb_rep", i, ".csv", sep = "")
  
  sim.1000kb <- read.table(path, sep = ",", header = T)
  
  sim.1000kb$thetaC <- sim.1000kb$MutRate - mean(sim.1000kb$MutRate)
  sim.1000kb$tmrcaC <- sim.1000kb$AverageTmrca - mean(sim.1000kb$AverageTmrca)
  m.div.1000kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.1000kb)
  anova.diversity <- Anova(m.div.1000kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.regular.1000kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.regular.1000kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.regular.1000kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.regular.1000kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

r2.regular.1r <- rbind.data.frame(t(colMeans(r2.regular.50kb)), t(colMeans(r2.regular.200kb)), t(colMeans(r2.regular.1000kb)))

# true landscapes, mu/r = 1

## 50 kb scale

r2.regular.50kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.regular.50kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.regular.50kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/RegularRegion/10r/RealLandscapes/RealLandscapes_50kb_rep", i, ".csv", sep = "")
  
  sim.50kb <- read.table(path, sep = ",", header = T)
  
  sim.50kb$thetaC <- sim.50kb$MutRate - mean(sim.50kb$MutRate)
  sim.50kb$tmrcaC <- sim.50kb$AverageTmrca - mean(sim.50kb$AverageTmrca)
  m.div.50kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.50kb)
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.regular.50kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.regular.50kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.regular.50kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.regular.50kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

## 200 kb scale

r2.regular.200kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.regular.200kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.regular.200kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/RegularRegion/10r/RealLandscapes/RealLandscapes_200kb_rep", i, ".csv", sep = "")
  
  sim.200kb <- read.table(path, sep = ",", header = T)
  
  sim.200kb$thetaC <- sim.200kb$MutRate - mean(sim.200kb$MutRate)
  sim.200kb$tmrcaC <- sim.200kb$AverageTmrca - mean(sim.200kb$AverageTmrca)
  m.div.200kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.200kb)
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.regular.200kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.regular.200kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.regular.200kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.regular.200kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

## 1000 kb scale

r2.regular.1000kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.regular.1000kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.regular.1000kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/RegularRegion/10r/RealLandscapes/RealLandscapes_1000kb_rep", i, ".csv", sep = "")
  
  sim.1000kb <- read.table(path, sep = ",", header = T)
  
  sim.1000kb$thetaC <- sim.1000kb$MutRate - mean(sim.1000kb$MutRate)
  sim.1000kb$tmrcaC <- sim.1000kb$AverageTmrca - mean(sim.1000kb$AverageTmrca)
  m.div.1000kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.1000kb)
  anova.diversity <- Anova(m.div.1000kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.regular.1000kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.regular.1000kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.regular.1000kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.regular.1000kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

r2.regular.10r <- rbind.data.frame(t(colMeans(r2.regular.50kb)), t(colMeans(r2.regular.200kb)), t(colMeans(r2.regular.1000kb)))

#Droso

# true landscapes, mu/r = 10

## 50 kb scale

r2.droso.50kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.droso.50kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.droso.50kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/Droso/1r/RealLandscapes/RealLandscapes_50kb_rep", i, ".csv", sep = "")
  
  sim.50kb <- read.table(path, sep = ",", header = T)
  
  sim.50kb$thetaC <- sim.50kb$MutRate - mean(sim.50kb$MutRate)
  sim.50kb$tmrcaC <- sim.50kb$AverageTmrca - mean(sim.50kb$AverageTmrca)
  m.div.50kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.50kb)
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.droso.50kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.droso.50kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.droso.50kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.droso.50kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

## 200 kb scale

r2.droso.200kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.droso.200kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.droso.200kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/Droso/1r/RealLandscapes/RealLandscapes_200kb_rep", i, ".csv", sep = "")
  
  sim.200kb <- read.table(path, sep = ",", header = T)
  
  sim.200kb$thetaC <- sim.200kb$MutRate - mean(sim.200kb$MutRate)
  sim.200kb$tmrcaC <- sim.200kb$AverageTmrca - mean(sim.200kb$AverageTmrca)
  m.div.200kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.200kb)
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.droso.200kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.droso.200kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.droso.200kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.droso.200kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

## 1000 kb scale

r2.droso.1000kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.droso.1000kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.droso.1000kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/Droso/1r/RealLandscapes/RealLandscapes_1000kb_rep", i, ".csv", sep = "")
  
  sim.1000kb <- read.table(path, sep = ",", header = T)
  
  sim.1000kb$thetaC <- sim.1000kb$MutRate - mean(sim.1000kb$MutRate)
  sim.1000kb$tmrcaC <- sim.1000kb$AverageTmrca - mean(sim.1000kb$AverageTmrca)
  m.div.1000kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.1000kb)
  anova.diversity <- Anova(m.div.1000kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.droso.1000kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.droso.1000kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.droso.1000kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.droso.1000kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

r2.droso.1r <- rbind.data.frame(t(colMeans(r2.droso.50kb)), t(colMeans(r2.droso.200kb)), t(colMeans(r2.droso.1000kb)))

# true landscapes, mu/r = 1

## 50 kb scale

r2.droso.50kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.droso.50kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.droso.50kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/Droso/10r/RealLandscapes/RealLandscapes_50kb_rep", i, ".csv", sep = "")
  
  sim.50kb <- read.table(path, sep = ",", header = T)
  
  sim.50kb$thetaC <- sim.50kb$MutRate - mean(sim.50kb$MutRate)
  sim.50kb$tmrcaC <- sim.50kb$AverageTmrca - mean(sim.50kb$AverageTmrca)
  m.div.50kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.50kb)
  anova.diversity <- Anova(m.div.50kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.droso.50kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.droso.50kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.droso.50kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.droso.50kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

## 200 kb scale

r2.droso.200kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.droso.200kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.droso.200kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/Droso/10r/RealLandscapes/RealLandscapes_200kb_rep", i, ".csv", sep = "")
  
  sim.200kb <- read.table(path, sep = ",", header = T)
  
  sim.200kb$thetaC <- sim.200kb$MutRate - mean(sim.200kb$MutRate)
  sim.200kb$tmrcaC <- sim.200kb$AverageTmrca - mean(sim.200kb$AverageTmrca)
  m.div.200kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.200kb)
  anova.diversity <- Anova(m.div.200kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.droso.200kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.droso.200kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.droso.200kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.droso.200kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

## 1000 kb scale

r2.droso.1000kb <- as.data.frame(matrix(nrow = nreps, ncol = 4))
colnames(r2.droso.1000kb) <- c("Total", "Theta", "TMRCA", "Theta:TMRCA")
row.names(r2.droso.1000kb) <- reps

for(i in reps) {
  
  path <- paste("bgs/Droso/10r/RealLandscapes/RealLandscapes_1000kb_rep", i, ".csv", sep = "")
  
  sim.1000kb <- read.table(path, sep = ",", header = T)
  
  sim.1000kb$thetaC <- sim.1000kb$MutRate - mean(sim.1000kb$MutRate)
  sim.1000kb$tmrcaC <- sim.1000kb$AverageTmrca - mean(sim.1000kb$AverageTmrca)
  m.div.1000kb <- lm(Pi.NonHomogeneous ~ thetaC + tmrcaC + thetaC:tmrcaC, data = sim.1000kb)
  anova.diversity <- Anova(m.div.1000kb)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  r2.droso.1000kb[i, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3]) * 100
  r2.droso.1000kb[i, 2] <- anova.diversity$VarExp[1] * 100
  r2.droso.1000kb[i, 3] <- anova.diversity$VarExp[2] * 100
  r2.droso.1000kb[i, 4] <- anova.diversity$VarExp[3] * 100
}

r2.droso.10r <- rbind.data.frame(t(colMeans(r2.droso.50kb)), t(colMeans(r2.droso.200kb)), t(colMeans(r2.droso.1000kb)))

r2.true.avg <- rbind.data.frame(r2.regular.1r, r2.regular.10r, r2.droso.1r, r2.droso.10r)
r2.true.avg$scale <- rep(c(50, 200, 1000), 4)
r2.true.avg$r <- as.factor(rep(c(rep("mu/r=10", 3), rep("mu/r=1", 3)), 2))
r2.true.avg$type <- c(rep("regular", 6), rep("droso", 6))

molten.r2 <- pivot_longer(r2.true.avg, cols = c("Total", "Theta", "TMRCA", "Theta:TMRCA"), names_to = "factor")
r2.plot <- ggplot(data = molten.r2, aes(x = scale, y = value, colour = factor))
r2.plot <- r2.plot + geom_line(data = molten.r2) + facet_grid(type~r)
r2.plot <- r2.plot + geom_point(aes(colour = factor), size = 4)
r2.plot <- r2.plot + scale_x_continuous(breaks = c(50, 200, 1000), trans="log10") 
r2.plot <- r2.plot + scale_y_continuous(breaks = pretty_breaks())
r2.plot <- r2.plot + labs(title = NULL, x = "Scale (kb)", y = "Variance Explained (%)") + theme_bw()
r2.plot <- r2.plot + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16))
r2.plot
ggsave("R2.bgs.pdf", r2.plot, device = "pdf", height = 12, width = 12)

