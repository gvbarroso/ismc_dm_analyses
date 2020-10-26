
library(reshape2)
library(tidyverse)
library(car)
library(scales)
library(cowplot)

################################################
#
# bottleneck_r_1e-8_flat_mu
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/bottleneck_r_1e-8_flat_mu/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.50kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)

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

sim.rho.200kb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_1$rho <- 1e-8
r2.sim.avg_1$mu_block <- "flat"

################################################
#
# bottleneck_r_1e-8_mu_change_50kb
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/bottleneck_r_1e-8_mu_change_50kb/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)
sim.theta.50kb <- read.table("sim.theta.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.2e+05.map", header = T)
sim.theta.200kb <- read.table("sim.theta.2e+05.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.1e+06.map", header = T)
sim.theta.1Mb <- read.table("sim.theta.1e+06.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_2$rho <- 1e-8
r2.sim.avg_2$mu_block <- "50 kb"


################################################
#
# bottleneck_r_1e-8_mu_change_500kb
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/bottleneck_r_1e-8_mu_change_500kb/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)
sim.theta.50kb <- read.table("sim.theta.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.2e+05.map", header = T)
sim.theta.200kb <- read.table("sim.theta.2e+05.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.1e+06.map", header = T)
sim.theta.1Mb <- read.table("sim.theta.1e+06.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_3$rho <- 1e-8
r2.sim.avg_3$mu_block <- "500 kb"


################################################
#
# bottleneck_r_1e-9_flat_mu
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/bottleneck_r_1e-9_flat_mu/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.50kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_4$rho <- 1e-9
r2.sim.avg_4$mu_block <- "flat"


################################################
#
# bottleneck_r_1e-9_mu_change_50kb
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/bottleneck_r_1e-9_mu_change_50kb/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)
sim.theta.50kb <- read.table("sim.theta.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.2e+05.map", header = T)
sim.theta.200kb <- read.table("sim.theta.2e+05.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.1e+06.map", header = T)
sim.theta.1Mb <- read.table("sim.theta.1e+06.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_5$rho <- 1e-9
r2.sim.avg_5$mu_block <- "50 kb"



################################################
#
# bottleneck_r_1e-9_mu_change_500kb
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/bottleneck_r_1e-9_mu_change_500kb/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)
sim.theta.50kb <- read.table("sim.theta.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.2e+05.map", header = T)
sim.theta.200kb <- read.table("sim.theta.2e+05.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.1e+06.map", header = T)
sim.theta.1Mb <- read.table("sim.theta.1e+06.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_6$rho <- 1e-9
r2.sim.avg_6$mu_block <- "500 kb"

################################################
#
# bottleneck plot
#
################################################

r2.sim.bottleneck <- rbind.data.frame(r2.sim.avg_1, r2.sim.avg_2, r2.sim.avg_3, r2.sim.avg_4, r2.sim.avg_5, r2.sim.avg_6, make.row.names = F)

molten.r2 <- melt(r2.sim.bottleneck, id.vars = c("bin.size", "rho", "mu_block"))

r2.plot_bottleneck <- ggplot(data = molten.r2, aes(x = bin.size, y = value, shape = variable))
r2.plot_bottleneck <- r2.plot_bottleneck + geom_line(data = molten.r2) + facet_grid(rho ~ mu_block) 
r2.plot_bottleneck <- r2.plot_bottleneck + geom_point(aes(shape = variable), size = 4)
r2.plot_bottleneck <- r2.plot_bottleneck + scale_shape_manual(values = c(0, 1, 2, 8))
r2.plot_bottleneck <- r2.plot_bottleneck + scale_x_continuous(breaks = c(50, 200, 1000), limits = c(40, 1200), trans="log10") 
r2.plot_bottleneck <- r2.plot_bottleneck + scale_y_continuous(breaks = pretty_breaks())
r2.plot_bottleneck <- r2.plot_bottleneck + labs(title = NULL, x = "Bin Size (kb)", y = "Variance Explained (%)") + theme_bw()
r2.plot_bottleneck <- r2.plot_bottleneck + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12))
ggsave("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/r2.bottleneck.pdf", r2.plot_bottleneck, device = "pdf")

################################################
#
# flat_Ne_r_1e-8_flat_mu
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/flat_Ne_r_1e-8_flat_mu/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.50kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_1$rho <- 1e-8
r2.sim.avg_1$mu_block <- "flat"

################################################
#
# flat_Ne_r_1e-8_mu_change_50kb
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/flat_Ne_r_1e-8_mu_change_50kb/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)
sim.theta.50kb <- read.table("sim.theta.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.2e+05.map", header = T)
sim.theta.200kb <- read.table("sim.theta.2e+05.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.1e+06.map", header = T)
sim.theta.1Mb <- read.table("sim.theta.1e+06.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_2$rho <- 1e-8
r2.sim.avg_2$mu_block <- "50 kb"


################################################
#
# flat_Ne_r_1e-8_mu_change_500kb
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/flat_Ne_r_1e-8_mu_change_500kb/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)
sim.theta.50kb <- read.table("sim.theta.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.2e+05.map", header = T)
sim.theta.200kb <- read.table("sim.theta.2e+05.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.1e+06.map", header = T)
sim.theta.1Mb <- read.table("sim.theta.1e+06.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_3$rho <- 1e-8
r2.sim.avg_3$mu_block <- "500 kb"


################################################
#
# flat_Ne_r_1e-9_flat_mu
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/flat_Ne_r_1e-9_flat_mu/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 3))
row.names(r2.sim.50kb) <- c("Total", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_4$rho <- 1e-9
r2.sim.avg_4$mu_block <- "flat"


################################################
#
# flat_Ne_r_1e-9_mu_change_50kb
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/flat_Ne_r_1e-9_mu_change_50kb/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)
sim.theta.50kb <- read.table("sim.theta.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.2e+05.map", header = T)
sim.theta.200kb <- read.table("sim.theta.2e+05.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.1e+06.map", header = T)
sim.theta.1Mb <- read.table("sim.theta.1e+06.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_5$rho <- 1e-9
r2.sim.avg_5$mu_block <- "50 kb"



################################################
#
# flat_Ne_r_1e-9_mu_change_500kb
#
################################################

setwd("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/flat_Ne_r_1e-9_mu_change_500kb/")

nreps <- 10
reps <- character(length = nreps)
for(i in 1:nreps) { reps[i] <- paste("rep_", i, sep = "") }

# 50kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 4))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA")
colnames(r2.sim.50kb) <- reps

sim.rho.50kb <- read.table("sim.rho.50000.map", header = T)
sim.theta.50kb <- read.table("sim.theta.50000.map", header = T)

for(i in 1:nreps) {
  rep_pi.50kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.50kb.bedgraph", sep = ""), header = T)
  rep_tmrca.50kb <- read.table(paste("rep_", i, "/sim.tmrca.50k.map", sep = ""), header = T)
  
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

sim.rho.200kb <- read.table("sim.rho.2e+05.map", header = T)
sim.theta.200kb <- read.table("sim.theta.2e+05.map", header = T)

for(i in 1:nreps) {
  rep_pi.200kb <- read.table(paste("rep_", i, "/rep_", i, ".pi.200kb.bedgraph", sep = ""), header = T)
  rep_tmrca.200kb <- read.table(paste("rep_", i, "/sim.tmrca.200k.map", sep = ""), header = T)
  
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

sim.rho.1Mb <- read.table("sim.rho.1e+06.map", header = T)
sim.theta.1Mb <- read.table("sim.theta.1e+06.map", header = T)

for(i in 1:nreps) {
  rep_pi.1Mb <- read.table(paste("rep_", i, "/rep_", i, ".pi.1Mb.bedgraph", sep = ""), header = T)
  rep_tmrca.1Mb <- read.table(paste("rep_", i, "/sim.tmrca.1M.map", sep = ""), header = T)
  
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
r2.sim.avg_6$rho <- 1e-9
r2.sim.avg_6$mu_block <- "500 kb"

################################################
#
# flat_Ne plot
#
################################################

r2.sim.flat_Ne <- rbind.data.frame(r2.sim.avg_1, r2.sim.avg_2, r2.sim.avg_3, r2.sim.avg_4, r2.sim.avg_5, r2.sim.avg_6, make.row.names = F)

molten.r2 <- melt(r2.sim.flat_Ne, id.vars = c("bin.size", "rho", "mu_block"))

r2.plot_flat_Ne <- ggplot(data = molten.r2, aes(x = bin.size, y = value, shape = variable))
r2.plot_flat_Ne <- r2.plot_flat_Ne + geom_line(data = molten.r2) + facet_grid(rho ~ mu_block) 
r2.plot_flat_Ne <- r2.plot_flat_Ne + geom_point(aes(shape = variable), size = 4)
r2.plot_flat_Ne <- r2.plot_flat_Ne + scale_shape_manual(values = c(0, 1, 2, 8))
r2.plot_flat_Ne <- r2.plot_flat_Ne + scale_x_continuous(breaks = c(50, 200, 1000), limits = c(40, 1200), trans="log10") 
r2.plot_flat_Ne <- r2.plot_flat_Ne + scale_y_continuous(breaks = pretty_breaks())
r2.plot_flat_Ne <- r2.plot_flat_Ne + labs(title = NULL, x = "Bin Size (kb)", y = "Variance Explained (%)") + theme_bw()
r2.plot_flat_Ne <- r2.plot_flat_Ne + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12))
ggsave("~/../../data2/gvbarroso/Data/iSMC/other_scenarios/r2.flat_Ne.pdf", r2.plot_flat_Ne, device = "pdf")
