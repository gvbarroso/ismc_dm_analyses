

library(MASS)
library(tidyverse)
library(car)
library(lmtest)
library(nlme)
library(cowplot)
library(scales)
library(magrittr)
library(reshape2)
library(scales)

setSessionTimeLimit(cpu = Inf, elapsed = Inf) # some of the plots in the script can take a few seconds to generate

########################
#
# var mu & true landscapes
#
########################

setwd("~/Data/iSMC/theta_paper/BGS_sims/var_mu/")

# 50 kb

R2.tab.50kb <- as.data.frame(matrix(nrow = 10, ncol = 5))
names(R2.tab.50kb) <- c("tmrca", "theta", "rho", "tmrca:theta", "total")

bgs.varmut.theta <- read.table("sim.theta.50000.map", header = T, sep = "\t")
bgs.varmut.rho <- read.table("sim.rho.50000.map", header = T, sep = "\t")

## rep 1
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_1_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_1.50kb.windowed.pi", header = T, sep = "\t")

rep1.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.50kb.maps$thetaC <- rep1.50kb.maps$theta - mean(rep1.50kb.maps$theta)
rep1.50kb.maps$tmrcaC <- rep1.50kb.maps$tmrca - mean(rep1.50kb.maps$tmrca)
rep1.50kb.maps$rhoC <- rep1.50kb.maps$rho - mean(rep1.50kb.maps$rho)

rep1.50kb.maps$bin <- 1:nrow(rep1.50kb.maps)

m.rep1.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.50kb.maps)

anova.diversity <- Anova(m.rep1.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[1,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 2
bgs.varmut.tmrca.rep2 <- read.table("TMRCA/dm2L_bgs_rep_2_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep2 <- read.table("Diversity/rep_2.50kb.windowed.pi", header = T, sep = "\t")

rep2.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep2$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep2$AverageTmrca)

names(rep2.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep2.50kb.maps$thetaC <- rep2.50kb.maps$theta - mean(rep2.50kb.maps$theta)
rep2.50kb.maps$tmrcaC <- rep2.50kb.maps$tmrca - mean(rep2.50kb.maps$tmrca)
rep2.50kb.maps$rhoC <- rep2.50kb.maps$rho - mean(rep2.50kb.maps$rho)

rep2.50kb.maps$bin <- 1:nrow(rep2.50kb.maps)

m.rep2.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep2.50kb.maps)

anova.diversity <- Anova(m.rep2.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[2,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 3
bgs.varmut.tmrca.rep3 <- read.table("TMRCA/dm2L_bgs_rep_3_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep3 <- read.table("Diversity/rep_3.50kb.windowed.pi", header = T, sep = "\t")

rep3.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep3$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep3$AverageTmrca)

names(rep3.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep3.50kb.maps$thetaC <- rep3.50kb.maps$theta - mean(rep3.50kb.maps$theta)
rep3.50kb.maps$tmrcaC <- rep3.50kb.maps$tmrca - mean(rep3.50kb.maps$tmrca)
rep3.50kb.maps$rhoC <- rep3.50kb.maps$rho - mean(rep3.50kb.maps$rho)

rep3.50kb.maps$bin <- 1:nrow(rep3.50kb.maps)

m.rep3.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep3.50kb.maps)

anova.diversity <- Anova(m.rep3.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[3,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 4
bgs.varmut.tmrca.rep4 <- read.table("TMRCA/dm2L_bgs_rep_4_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep4 <- read.table("Diversity/rep_4.50kb.windowed.pi", header = T, sep = "\t")

rep4.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep4$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep4$AverageTmrca)

names(rep4.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep4.50kb.maps$thetaC <- rep4.50kb.maps$theta - mean(rep4.50kb.maps$theta)
rep4.50kb.maps$tmrcaC <- rep4.50kb.maps$tmrca - mean(rep4.50kb.maps$tmrca)
rep4.50kb.maps$rhoC <- rep4.50kb.maps$rho - mean(rep4.50kb.maps$rho)

rep4.50kb.maps$bin <- 1:nrow(rep4.50kb.maps)

m.rep4.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep4.50kb.maps)

anova.diversity <- Anova(m.rep4.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[4,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 5
bgs.varmut.tmrca.rep5 <- read.table("TMRCA/dm2L_bgs_rep_5_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep5 <- read.table("Diversity/rep_5.50kb.windowed.pi", header = T, sep = "\t")

rep5.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep5$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep5$AverageTmrca)

names(rep5.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep5.50kb.maps$thetaC <- rep5.50kb.maps$theta - mean(rep5.50kb.maps$theta)
rep5.50kb.maps$tmrcaC <- rep5.50kb.maps$tmrca - mean(rep5.50kb.maps$tmrca)
rep5.50kb.maps$rhoC <- rep5.50kb.maps$rho - mean(rep5.50kb.maps$rho)

rep5.50kb.maps$bin <- 1:nrow(rep5.50kb.maps)

m.rep5.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep5.50kb.maps)

anova.diversity <- Anova(m.rep5.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[5,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 6
bgs.varmut.tmrca.rep6 <- read.table("TMRCA/dm2L_bgs_rep_6_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep6 <- read.table("Diversity/rep_6.50kb.windowed.pi", header = T, sep = "\t")

rep6.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep6$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep6$AverageTmrca)

names(rep6.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep6.50kb.maps$thetaC <- rep6.50kb.maps$theta - mean(rep6.50kb.maps$theta)
rep6.50kb.maps$tmrcaC <- rep6.50kb.maps$tmrca - mean(rep6.50kb.maps$tmrca)
rep6.50kb.maps$rhoC <- rep6.50kb.maps$rho - mean(rep6.50kb.maps$rho)

rep6.50kb.maps$bin <- 1:nrow(rep6.50kb.maps)

m.rep6.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep6.50kb.maps)

anova.diversity <- Anova(m.rep6.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[6,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 7
bgs.varmut.tmrca.rep7 <- read.table("TMRCA/dm2L_bgs_rep_7_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep7 <- read.table("Diversity/rep_7.50kb.windowed.pi", header = T, sep = "\t")

rep7.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep7$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep7$AverageTmrca)

names(rep7.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep7.50kb.maps$thetaC <- rep7.50kb.maps$theta - mean(rep7.50kb.maps$theta)
rep7.50kb.maps$tmrcaC <- rep7.50kb.maps$tmrca - mean(rep7.50kb.maps$tmrca)
rep7.50kb.maps$rhoC <- rep7.50kb.maps$rho - mean(rep7.50kb.maps$rho)

rep7.50kb.maps$bin <- 1:nrow(rep7.50kb.maps)

m.rep7.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep7.50kb.maps)

anova.diversity <- Anova(m.rep7.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[7,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 8
bgs.varmut.tmrca.rep8 <- read.table("TMRCA/dm2L_bgs_rep_8_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep8 <- read.table("Diversity/rep_8.50kb.windowed.pi", header = T, sep = "\t")

rep8.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep8$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep8$AverageTmrca)

names(rep8.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep8.50kb.maps$thetaC <- rep8.50kb.maps$theta - mean(rep8.50kb.maps$theta)
rep8.50kb.maps$tmrcaC <- rep8.50kb.maps$tmrca - mean(rep8.50kb.maps$tmrca)
rep8.50kb.maps$rhoC <- rep8.50kb.maps$rho - mean(rep8.50kb.maps$rho)

rep8.50kb.maps$bin <- 1:nrow(rep8.50kb.maps)

m.rep8.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep8.50kb.maps)

anova.diversity <- Anova(m.rep8.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[8,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 9
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_9_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_9.50kb.windowed.pi", header = T, sep = "\t")

rep1.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.50kb.maps$thetaC <- rep1.50kb.maps$theta - mean(rep1.50kb.maps$theta)
rep1.50kb.maps$tmrcaC <- rep1.50kb.maps$tmrca - mean(rep1.50kb.maps$tmrca)
rep1.50kb.maps$rhoC <- rep1.50kb.maps$rho - mean(rep1.50kb.maps$rho)

rep1.50kb.maps$bin <- 1:nrow(rep1.50kb.maps)

m.rep1.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.50kb.maps)

anova.diversity <- Anova(m.rep1.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[9,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)


## rep 10
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_10_tmrca_50kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_10.50kb.windowed.pi", header = T, sep = "\t")

rep1.50kb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.50kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.50kb.maps$thetaC <- rep1.50kb.maps$theta - mean(rep1.50kb.maps$theta)
rep1.50kb.maps$tmrcaC <- rep1.50kb.maps$tmrca - mean(rep1.50kb.maps$tmrca)
rep1.50kb.maps$rhoC <- rep1.50kb.maps$rho - mean(rep1.50kb.maps$rho)

rep1.50kb.maps$bin <- 1:nrow(rep1.50kb.maps)

m.rep1.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.50kb.maps)

anova.diversity <- Anova(m.rep1.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.50kb[10,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)


# 200 kb

R2.tab.200kb <- as.data.frame(matrix(nrow = 10, ncol = 5))
names(R2.tab.200kb) <- c("tmrca", "theta", "rho", "tmrca:theta", "total")

bgs.varmut.theta <- read.table("sim.theta.2e+05.map", header = T, sep = "\t")
bgs.varmut.rho <- read.table("sim.rho.2e+05.map", header = T, sep = "\t")

## rep 1
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_1_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_1.200kb.windowed.pi", header = T, sep = "\t")

rep1.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.200kb.maps$thetaC <- rep1.200kb.maps$theta - mean(rep1.200kb.maps$theta)
rep1.200kb.maps$tmrcaC <- rep1.200kb.maps$tmrca - mean(rep1.200kb.maps$tmrca)
rep1.200kb.maps$rhoC <- rep1.200kb.maps$rho - mean(rep1.200kb.maps$rho)

rep1.200kb.maps$bin <- 1:nrow(rep1.200kb.maps)

m.rep1.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.200kb.maps)

anova.diversity <- Anova(m.rep1.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[1,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 2
bgs.varmut.tmrca.rep2 <- read.table("TMRCA/dm2L_bgs_rep_2_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep2 <- read.table("Diversity/rep_2.200kb.windowed.pi", header = T, sep = "\t")

rep2.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep2$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep2$AverageTmrca)

names(rep2.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep2.200kb.maps$thetaC <- rep2.200kb.maps$theta - mean(rep2.200kb.maps$theta)
rep2.200kb.maps$tmrcaC <- rep2.200kb.maps$tmrca - mean(rep2.200kb.maps$tmrca)
rep2.200kb.maps$rhoC <- rep2.200kb.maps$rho - mean(rep2.200kb.maps$rho)

rep2.200kb.maps$bin <- 1:nrow(rep2.200kb.maps)

m.rep2.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep2.200kb.maps)

anova.diversity <- Anova(m.rep2.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[2,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 3
bgs.varmut.tmrca.rep3 <- read.table("TMRCA/dm2L_bgs_rep_3_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep3 <- read.table("Diversity/rep_3.200kb.windowed.pi", header = T, sep = "\t")

rep3.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep3$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep3$AverageTmrca)

names(rep3.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep3.200kb.maps$thetaC <- rep3.200kb.maps$theta - mean(rep3.200kb.maps$theta)
rep3.200kb.maps$tmrcaC <- rep3.200kb.maps$tmrca - mean(rep3.200kb.maps$tmrca)
rep3.200kb.maps$rhoC <- rep3.200kb.maps$rho - mean(rep3.200kb.maps$rho)

rep3.200kb.maps$bin <- 1:nrow(rep3.200kb.maps)

m.rep3.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep3.200kb.maps)

anova.diversity <- Anova(m.rep3.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[3,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 4
bgs.varmut.tmrca.rep4 <- read.table("TMRCA/dm2L_bgs_rep_4_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep4 <- read.table("Diversity/rep_4.200kb.windowed.pi", header = T, sep = "\t")

rep4.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep4$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep4$AverageTmrca)

names(rep4.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep4.200kb.maps$thetaC <- rep4.200kb.maps$theta - mean(rep4.200kb.maps$theta)
rep4.200kb.maps$tmrcaC <- rep4.200kb.maps$tmrca - mean(rep4.200kb.maps$tmrca)
rep4.200kb.maps$rhoC <- rep4.200kb.maps$rho - mean(rep4.200kb.maps$rho)

rep4.200kb.maps$bin <- 1:nrow(rep4.200kb.maps)

m.rep4.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep4.200kb.maps)

anova.diversity <- Anova(m.rep4.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[4,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 5
bgs.varmut.tmrca.rep5 <- read.table("TMRCA/dm2L_bgs_rep_5_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep5 <- read.table("Diversity/rep_5.200kb.windowed.pi", header = T, sep = "\t")

rep5.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep5$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep5$AverageTmrca)

names(rep5.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep5.200kb.maps$thetaC <- rep5.200kb.maps$theta - mean(rep5.200kb.maps$theta)
rep5.200kb.maps$tmrcaC <- rep5.200kb.maps$tmrca - mean(rep5.200kb.maps$tmrca)
rep5.200kb.maps$rhoC <- rep5.200kb.maps$rho - mean(rep5.200kb.maps$rho)

rep5.200kb.maps$bin <- 1:nrow(rep5.200kb.maps)

m.rep5.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep5.200kb.maps)

anova.diversity <- Anova(m.rep5.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[5,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 6
bgs.varmut.tmrca.rep6 <- read.table("TMRCA/dm2L_bgs_rep_6_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep6 <- read.table("Diversity/rep_6.200kb.windowed.pi", header = T, sep = "\t")

rep6.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep6$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep6$AverageTmrca)

names(rep6.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep6.200kb.maps$thetaC <- rep6.200kb.maps$theta - mean(rep6.200kb.maps$theta)
rep6.200kb.maps$tmrcaC <- rep6.200kb.maps$tmrca - mean(rep6.200kb.maps$tmrca)
rep6.200kb.maps$rhoC <- rep6.200kb.maps$rho - mean(rep6.200kb.maps$rho)

rep6.200kb.maps$bin <- 1:nrow(rep6.200kb.maps)

m.rep6.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep6.200kb.maps)

anova.diversity <- Anova(m.rep6.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[6,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 7
bgs.varmut.tmrca.rep7 <- read.table("TMRCA/dm2L_bgs_rep_7_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep7 <- read.table("Diversity/rep_7.200kb.windowed.pi", header = T, sep = "\t")

rep7.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep7$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep7$AverageTmrca)

names(rep7.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep7.200kb.maps$thetaC <- rep7.200kb.maps$theta - mean(rep7.200kb.maps$theta)
rep7.200kb.maps$tmrcaC <- rep7.200kb.maps$tmrca - mean(rep7.200kb.maps$tmrca)
rep7.200kb.maps$rhoC <- rep7.200kb.maps$rho - mean(rep7.200kb.maps$rho)

rep7.200kb.maps$bin <- 1:nrow(rep7.200kb.maps)

m.rep7.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep7.200kb.maps)

anova.diversity <- Anova(m.rep7.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[7,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 8
bgs.varmut.tmrca.rep8 <- read.table("TMRCA/dm2L_bgs_rep_8_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep8 <- read.table("Diversity/rep_8.200kb.windowed.pi", header = T, sep = "\t")

rep8.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep8$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep8$AverageTmrca)

names(rep8.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep8.200kb.maps$thetaC <- rep8.200kb.maps$theta - mean(rep8.200kb.maps$theta)
rep8.200kb.maps$tmrcaC <- rep8.200kb.maps$tmrca - mean(rep8.200kb.maps$tmrca)
rep8.200kb.maps$rhoC <- rep8.200kb.maps$rho - mean(rep8.200kb.maps$rho)

rep8.200kb.maps$bin <- 1:nrow(rep8.200kb.maps)

m.rep8.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep8.200kb.maps)

anova.diversity <- Anova(m.rep8.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[8,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)

## rep 9
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_9_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_9.200kb.windowed.pi", header = T, sep = "\t")

rep1.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.200kb.maps$thetaC <- rep1.200kb.maps$theta - mean(rep1.200kb.maps$theta)
rep1.200kb.maps$tmrcaC <- rep1.200kb.maps$tmrca - mean(rep1.200kb.maps$tmrca)
rep1.200kb.maps$rhoC <- rep1.200kb.maps$rho - mean(rep1.200kb.maps$rho)

rep1.200kb.maps$bin <- 1:nrow(rep1.200kb.maps)

m.rep1.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.200kb.maps)

anova.diversity <- Anova(m.rep1.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[9,] <- c(anova.diversity$VarExp[3] * 100, 
                     anova.diversity$VarExp[1] * 100,
                     anova.diversity$VarExp[2] * 100,
                     anova.diversity$VarExp[4] * 100,
                     100 - anova.diversity$VarExp[5] * 100)


## rep 10
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_10_tmrca_200kb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_10.200kb.windowed.pi", header = T, sep = "\t")

rep1.200kb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                   bgs.varmut.theta$sim,
                                   bgs.varmut.rho$sim,
                                   bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.200kb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.200kb.maps$thetaC <- rep1.200kb.maps$theta - mean(rep1.200kb.maps$theta)
rep1.200kb.maps$tmrcaC <- rep1.200kb.maps$tmrca - mean(rep1.200kb.maps$tmrca)
rep1.200kb.maps$rhoC <- rep1.200kb.maps$rho - mean(rep1.200kb.maps$rho)

rep1.200kb.maps$bin <- 1:nrow(rep1.200kb.maps)

m.rep1.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.200kb.maps)

anova.diversity <- Anova(m.rep1.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.200kb[10,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

# 1 Mb

R2.tab.1Mb <- as.data.frame(matrix(nrow = 10, ncol = 5))
names(R2.tab.1Mb) <- c("tmrca", "theta", "rho", "tmrca:theta", "total")

bgs.varmut.theta <- read.table("sim.theta.1e+06.map", header = T, sep = "\t")
bgs.varmut.rho <- read.table("sim.rho.1e+06.map", header = T, sep = "\t")

## rep 1
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_1_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_1.1Mb.windowed.pi", header = T, sep = "\t")

rep1.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.1Mb.maps$thetaC <- rep1.1Mb.maps$theta - mean(rep1.1Mb.maps$theta)
rep1.1Mb.maps$tmrcaC <- rep1.1Mb.maps$tmrca - mean(rep1.1Mb.maps$tmrca)
rep1.1Mb.maps$rhoC <- rep1.1Mb.maps$rho - mean(rep1.1Mb.maps$rho)

rep1.1Mb.maps$bin <- 1:nrow(rep1.1Mb.maps)

m.rep1.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.1Mb.maps)

anova.diversity <- Anova(m.rep1.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[1,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

## rep 2
bgs.varmut.tmrca.rep2 <- read.table("TMRCA/dm2L_bgs_rep_2_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep2 <- read.table("Diversity/rep_2.1Mb.windowed.pi", header = T, sep = "\t")

rep2.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep2$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep2$AverageTmrca)

names(rep2.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep2.1Mb.maps$thetaC <- rep2.1Mb.maps$theta - mean(rep2.1Mb.maps$theta)
rep2.1Mb.maps$tmrcaC <- rep2.1Mb.maps$tmrca - mean(rep2.1Mb.maps$tmrca)
rep2.1Mb.maps$rhoC <- rep2.1Mb.maps$rho - mean(rep2.1Mb.maps$rho)

rep2.1Mb.maps$bin <- 1:nrow(rep2.1Mb.maps)

m.rep2.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep2.1Mb.maps)

anova.diversity <- Anova(m.rep2.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[2,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

## rep 3
bgs.varmut.tmrca.rep3 <- read.table("TMRCA/dm2L_bgs_rep_3_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep3 <- read.table("Diversity/rep_3.1Mb.windowed.pi", header = T, sep = "\t")

rep3.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep3$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep3$AverageTmrca)

names(rep3.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep3.1Mb.maps$thetaC <- rep3.1Mb.maps$theta - mean(rep3.1Mb.maps$theta)
rep3.1Mb.maps$tmrcaC <- rep3.1Mb.maps$tmrca - mean(rep3.1Mb.maps$tmrca)
rep3.1Mb.maps$rhoC <- rep3.1Mb.maps$rho - mean(rep3.1Mb.maps$rho)

rep3.1Mb.maps$bin <- 1:nrow(rep3.1Mb.maps)

m.rep3.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep3.1Mb.maps)

anova.diversity <- Anova(m.rep3.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[3,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

## rep 4
bgs.varmut.tmrca.rep4 <- read.table("TMRCA/dm2L_bgs_rep_4_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep4 <- read.table("Diversity/rep_4.1Mb.windowed.pi", header = T, sep = "\t")

rep4.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep4$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep4$AverageTmrca)

names(rep4.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep4.1Mb.maps$thetaC <- rep4.1Mb.maps$theta - mean(rep4.1Mb.maps$theta)
rep4.1Mb.maps$tmrcaC <- rep4.1Mb.maps$tmrca - mean(rep4.1Mb.maps$tmrca)
rep4.1Mb.maps$rhoC <- rep4.1Mb.maps$rho - mean(rep4.1Mb.maps$rho)

rep4.1Mb.maps$bin <- 1:nrow(rep4.1Mb.maps)

m.rep4.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep4.1Mb.maps)

anova.diversity <- Anova(m.rep4.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[4,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

## rep 5
bgs.varmut.tmrca.rep5 <- read.table("TMRCA/dm2L_bgs_rep_5_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep5 <- read.table("Diversity/rep_5.1Mb.windowed.pi", header = T, sep = "\t")

rep5.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep5$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep5$AverageTmrca)

names(rep5.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep5.1Mb.maps$thetaC <- rep5.1Mb.maps$theta - mean(rep5.1Mb.maps$theta)
rep5.1Mb.maps$tmrcaC <- rep5.1Mb.maps$tmrca - mean(rep5.1Mb.maps$tmrca)
rep5.1Mb.maps$rhoC <- rep5.1Mb.maps$rho - mean(rep5.1Mb.maps$rho)

rep5.1Mb.maps$bin <- 1:nrow(rep5.1Mb.maps)

m.rep5.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep5.1Mb.maps)

anova.diversity <- Anova(m.rep5.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[5,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

## rep 6
bgs.varmut.tmrca.rep6 <- read.table("TMRCA/dm2L_bgs_rep_6_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep6 <- read.table("Diversity/rep_6.1Mb.windowed.pi", header = T, sep = "\t")

rep6.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep6$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep6$AverageTmrca)

names(rep6.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep6.1Mb.maps$thetaC <- rep6.1Mb.maps$theta - mean(rep6.1Mb.maps$theta)
rep6.1Mb.maps$tmrcaC <- rep6.1Mb.maps$tmrca - mean(rep6.1Mb.maps$tmrca)
rep6.1Mb.maps$rhoC <- rep6.1Mb.maps$rho - mean(rep6.1Mb.maps$rho)

rep6.1Mb.maps$bin <- 1:nrow(rep6.1Mb.maps)

m.rep6.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep6.1Mb.maps)

anova.diversity <- Anova(m.rep6.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[6,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

## rep 7
bgs.varmut.tmrca.rep7 <- read.table("TMRCA/dm2L_bgs_rep_7_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep7 <- read.table("Diversity/rep_7.1Mb.windowed.pi", header = T, sep = "\t")

rep7.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep7$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep7$AverageTmrca)

names(rep7.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep7.1Mb.maps$thetaC <- rep7.1Mb.maps$theta - mean(rep7.1Mb.maps$theta)
rep7.1Mb.maps$tmrcaC <- rep7.1Mb.maps$tmrca - mean(rep7.1Mb.maps$tmrca)
rep7.1Mb.maps$rhoC <- rep7.1Mb.maps$rho - mean(rep7.1Mb.maps$rho)

rep7.1Mb.maps$bin <- 1:nrow(rep7.1Mb.maps)

m.rep7.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep7.1Mb.maps)

anova.diversity <- Anova(m.rep7.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[7,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

## rep 8
bgs.varmut.tmrca.rep8 <- read.table("TMRCA/dm2L_bgs_rep_8_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep8 <- read.table("Diversity/rep_8.1Mb.windowed.pi", header = T, sep = "\t")

rep8.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep8$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep8$AverageTmrca)

names(rep8.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep8.1Mb.maps$thetaC <- rep8.1Mb.maps$theta - mean(rep8.1Mb.maps$theta)
rep8.1Mb.maps$tmrcaC <- rep8.1Mb.maps$tmrca - mean(rep8.1Mb.maps$tmrca)
rep8.1Mb.maps$rhoC <- rep8.1Mb.maps$rho - mean(rep8.1Mb.maps$rho)

rep8.1Mb.maps$bin <- 1:nrow(rep8.1Mb.maps)

m.rep8.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep8.1Mb.maps)

anova.diversity <- Anova(m.rep8.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[8,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)

## rep 9
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_9_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_9.1Mb.windowed.pi", header = T, sep = "\t")

rep1.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.1Mb.maps$thetaC <- rep1.1Mb.maps$theta - mean(rep1.1Mb.maps$theta)
rep1.1Mb.maps$tmrcaC <- rep1.1Mb.maps$tmrca - mean(rep1.1Mb.maps$tmrca)
rep1.1Mb.maps$rhoC <- rep1.1Mb.maps$rho - mean(rep1.1Mb.maps$rho)

rep1.1Mb.maps$bin <- 1:nrow(rep1.1Mb.maps)

m.rep1.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.1Mb.maps)

anova.diversity <- Anova(m.rep1.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[9,] <- c(anova.diversity$VarExp[3] * 100, 
                      anova.diversity$VarExp[1] * 100,
                      anova.diversity$VarExp[2] * 100,
                      anova.diversity$VarExp[4] * 100,
                      100 - anova.diversity$VarExp[5] * 100)


## rep 10
bgs.varmut.tmrca.rep1 <- read.table("TMRCA/dm2L_bgs_rep_10_tmrca_1Mb.csv", header = T, sep = ",")
bgs.varmut.pi.rep1 <- read.table("Diversity/rep_10.1Mb.windowed.pi", header = T, sep = "\t")

rep1.1Mb.maps <- cbind.data.frame(bgs.varmut.pi.rep1$PI,
                                    bgs.varmut.theta$sim,
                                    bgs.varmut.rho$sim,
                                    bgs.varmut.tmrca.rep1$AverageTmrca)

names(rep1.1Mb.maps) <- c("diversity", "theta", "rho", "tmrca")

# centering
rep1.1Mb.maps$thetaC <- rep1.1Mb.maps$theta - mean(rep1.1Mb.maps$theta)
rep1.1Mb.maps$tmrcaC <- rep1.1Mb.maps$tmrca - mean(rep1.1Mb.maps$tmrca)
rep1.1Mb.maps$rhoC <- rep1.1Mb.maps$rho - mean(rep1.1Mb.maps$rho)

rep1.1Mb.maps$bin <- 1:nrow(rep1.1Mb.maps)

m.rep1.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep1.1Mb.maps)

anova.diversity <- Anova(m.rep1.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

R2.tab.1Mb[10,] <- c(anova.diversity$VarExp[3] * 100, 
                       anova.diversity$VarExp[1] * 100,
                       anova.diversity$VarExp[2] * 100,
                       anova.diversity$VarExp[4] * 100,
                       100 - anova.diversity$VarExp[5] * 100)

R2.tab <- rbind.data.frame(R2.tab.50kb, R2.tab.200kb, R2.tab.1Mb)
R2.tab$bin.size <- c(rep(5e+4, 10), rep(2e+5, 10), rep(1e+6, 10))

#########################
#
# flat mu & inferred landscapes
#
#########################

setwd("~/Data/iSMC/theta_paper/BGS_sims/flat_mu/")

# for plotting exons later on:
dm_2L_exome <- read.table("~/Data/iSMC/theta_paper/slim_sims/dm_tbl.txt", header = T) %>% 
  filter(seqnames == "1") %>%
  filter(!is.na(ensembl_exon_id)) %>%
  dplyr::select(start, end)

names(dm_2L_exome) <- c("exon_start", "exon_end")

# 50kb

## rep_1

rep_1.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_1.diversity.50kb.bedgraph", header = T)
rep_1.pi.50kb$sample_mean <- apply(rep_1.pi.50kb[,(4:ncol(rep_1.pi.50kb))], 1, mean)

rep_1.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_1.TMRCA.50kb.bedgraph", header = T)
rep_1.tmrca.50kb$sample_mean <- apply(rep_1.tmrca.50kb[4:ncol(rep_1.tmrca.50kb)], 1, mean)

rep_1.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_1.theta.50kb.bedgraph", header = T)
rep_1.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_1.rho.50kb.bedgraph", header = T)

rep_1.maps.50kb <- as.data.frame(cbind(rep_1.pi.50kb$sample_mean,
                                       rep_1.tmrca.50kb$sample_mean,
                                       rep_1.theta.50kb$sample_mean,
                                       rep_1.rho.50kb$sample_mean))

names(rep_1.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_1.maps.50kb$thetaC <- rep_1.maps.50kb$theta - mean(rep_1.maps.50kb$theta)
rep_1.maps.50kb$tmrcaC <- rep_1.maps.50kb$tmrca - mean(rep_1.maps.50kb$tmrca)
rep_1.maps.50kb$rhoC <- rep_1.maps.50kb$rho - mean(rep_1.maps.50kb$rho)

rep_1.maps.50kb$bin <- 1:nrow(rep_1.maps.50kb)

m.rep_1.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_1.maps.50kb)

rep_1.anova <- Anova(m.rep_1.50kb)
apiss <- rep_1.anova$"Sum Sq"
rep_1.anova$VarExp <- apiss / sum(apiss)



## rep_2

rep_2.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_2.diversity.50kb.bedgraph", header = T)
rep_2.pi.50kb$sample_mean <- apply(rep_2.pi.50kb[,(4:ncol(rep_2.pi.50kb))], 1, mean)

rep_2.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_2.TMRCA.50kb.bedgraph", header = T)
rep_2.tmrca.50kb$sample_mean <- apply(rep_2.tmrca.50kb[4:ncol(rep_2.tmrca.50kb)], 1, mean)

rep_2.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_2.theta.50kb.bedgraph", header = T)
rep_2.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_2.rho.50kb.bedgraph", header = T)

rep_2.maps.50kb <- as.data.frame(cbind(rep_2.pi.50kb$sample_mean,
                                       rep_2.tmrca.50kb$sample_mean,
                                       rep_2.theta.50kb$sample_mean,
                                       rep_2.rho.50kb$sample_mean))

names(rep_2.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_2.maps.50kb$thetaC <- rep_2.maps.50kb$theta - mean(rep_2.maps.50kb$theta)
rep_2.maps.50kb$tmrcaC <- rep_2.maps.50kb$tmrca - mean(rep_2.maps.50kb$tmrca)
rep_2.maps.50kb$rhoC <- rep_2.maps.50kb$rho - mean(rep_2.maps.50kb$rho)

rep_2.maps.50kb$bin <- 1:nrow(rep_2.maps.50kb)

m.rep_2.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_2.maps.50kb)

rep_2.anova <- Anova(m.rep_2.50kb)
apiss <- rep_2.anova$"Sum Sq"
rep_2.anova$VarExp <- apiss / sum(apiss)


## rep_3

rep_3.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_3.diversity.50kb.bedgraph", header = T)
rep_3.pi.50kb$sample_mean <- apply(rep_3.pi.50kb[,(4:ncol(rep_3.pi.50kb))], 1, mean)

rep_3.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_3.TMRCA.50kb.bedgraph", header = T)
rep_3.tmrca.50kb$sample_mean <- apply(rep_3.tmrca.50kb[4:ncol(rep_3.tmrca.50kb)], 1, mean)

rep_3.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_3.theta.50kb.bedgraph", header = T)
rep_3.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_3.rho.50kb.bedgraph", header = T)

rep_3.maps.50kb <- as.data.frame(cbind(rep_3.pi.50kb$sample_mean,
                                       rep_3.tmrca.50kb$sample_mean,
                                       rep_3.theta.50kb$sample_mean,
                                       rep_3.rho.50kb$sample_mean))

names(rep_3.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_3.maps.50kb$thetaC <- rep_3.maps.50kb$theta - mean(rep_3.maps.50kb$theta)
rep_3.maps.50kb$tmrcaC <- rep_3.maps.50kb$tmrca - mean(rep_3.maps.50kb$tmrca)
rep_3.maps.50kb$rhoC <- rep_3.maps.50kb$rho - mean(rep_3.maps.50kb$rho)

rep_3.maps.50kb$bin <- 1:nrow(rep_3.maps.50kb)

m.rep_3.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_3.maps.50kb)

rep_3.anova <- Anova(m.rep_3.50kb)
apiss <- rep_3.anova$"Sum Sq"
rep_3.anova$VarExp <- apiss / sum(apiss)



## rep_4

rep_4.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_4.diversity.50kb.bedgraph", header = T)
rep_4.pi.50kb$sample_mean <- apply(rep_4.pi.50kb[,(4:ncol(rep_4.pi.50kb))], 1, mean)

rep_4.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_4.TMRCA.50kb.bedgraph", header = T)
rep_4.tmrca.50kb$sample_mean <- apply(rep_4.tmrca.50kb[4:ncol(rep_4.tmrca.50kb)], 1, mean)

rep_4.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_4.theta.50kb.bedgraph", header = T)
rep_4.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_4.rho.50kb.bedgraph", header = T)

rep_4.maps.50kb <- as.data.frame(cbind(rep_4.pi.50kb$sample_mean,
                                       rep_4.tmrca.50kb$sample_mean,
                                       rep_4.theta.50kb$sample_mean,
                                       rep_4.rho.50kb$sample_mean))

names(rep_4.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_4.maps.50kb$thetaC <- rep_4.maps.50kb$theta - mean(rep_4.maps.50kb$theta)
rep_4.maps.50kb$tmrcaC <- rep_4.maps.50kb$tmrca - mean(rep_4.maps.50kb$tmrca)
rep_4.maps.50kb$rhoC <- rep_4.maps.50kb$rho - mean(rep_4.maps.50kb$rho)

rep_4.maps.50kb$bin <- 1:nrow(rep_4.maps.50kb)

m.rep_4.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_4.maps.50kb)

rep_4.anova <- Anova(m.rep_4.50kb)
apiss <- rep_4.anova$"Sum Sq"
rep_4.anova$VarExp <- apiss / sum(apiss)



## rep_5

rep_5.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_5.diversity.50kb.bedgraph", header = T)
rep_5.pi.50kb$sample_mean <- apply(rep_5.pi.50kb[,(4:ncol(rep_5.pi.50kb))], 1, mean)

rep_5.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_5.TMRCA.50kb.bedgraph", header = T)
rep_5.tmrca.50kb$sample_mean <- apply(rep_5.tmrca.50kb[4:ncol(rep_5.tmrca.50kb)], 1, mean)

rep_5.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_5.theta.50kb.bedgraph", header = T)
rep_5.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_5.rho.50kb.bedgraph", header = T)

rep_5.maps.50kb <- as.data.frame(cbind(rep_5.pi.50kb$sample_mean,
                                       rep_5.tmrca.50kb$sample_mean,
                                       rep_5.theta.50kb$sample_mean,
                                       rep_5.rho.50kb$sample_mean))

names(rep_5.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_5.maps.50kb$thetaC <- rep_5.maps.50kb$theta - mean(rep_5.maps.50kb$theta)
rep_5.maps.50kb$tmrcaC <- rep_5.maps.50kb$tmrca - mean(rep_5.maps.50kb$tmrca)
rep_5.maps.50kb$rhoC <- rep_5.maps.50kb$rho - mean(rep_5.maps.50kb$rho)

rep_5.maps.50kb$bin <- 1:nrow(rep_5.maps.50kb)

m.rep_5.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_5.maps.50kb)

rep_5.anova <- Anova(m.rep_5.50kb)
apiss <- rep_5.anova$"Sum Sq"
rep_5.anova$VarExp <- apiss / sum(apiss)



## rep_6

rep_6.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_6.diversity.50kb.bedgraph", header = T)
rep_6.pi.50kb$sample_mean <- apply(rep_6.pi.50kb[,(4:ncol(rep_6.pi.50kb))], 1, mean)

rep_6.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_6.TMRCA.50kb.bedgraph", header = T)
rep_6.tmrca.50kb$sample_mean <- apply(rep_6.tmrca.50kb[4:ncol(rep_6.tmrca.50kb)], 1, mean)

rep_6.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_6.theta.50kb.bedgraph", header = T)
rep_6.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_6.rho.50kb.bedgraph", header = T)

rep_6.maps.50kb <- as.data.frame(cbind(rep_6.pi.50kb$sample_mean,
                                       rep_6.tmrca.50kb$sample_mean,
                                       rep_6.theta.50kb$sample_mean,
                                       rep_6.rho.50kb$sample_mean))

names(rep_6.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_6.maps.50kb$thetaC <- rep_6.maps.50kb$theta - mean(rep_6.maps.50kb$theta)
rep_6.maps.50kb$tmrcaC <- rep_6.maps.50kb$tmrca - mean(rep_6.maps.50kb$tmrca)
rep_6.maps.50kb$rhoC <- rep_6.maps.50kb$rho - mean(rep_6.maps.50kb$rho)

rep_6.maps.50kb$bin <- 1:nrow(rep_6.maps.50kb)

m.rep_6.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_6.maps.50kb)

rep_6.anova <- Anova(m.rep_6.50kb)
apiss <- rep_6.anova$"Sum Sq"
rep_6.anova$VarExp <- apiss / sum(apiss)




## rep_7

rep_7.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_7.diversity.50kb.bedgraph", header = T)
rep_7.pi.50kb$sample_mean <- apply(rep_7.pi.50kb[,(4:ncol(rep_7.pi.50kb))], 1, mean)

rep_7.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_7.TMRCA.50kb.bedgraph", header = T)
rep_7.tmrca.50kb$sample_mean <- apply(rep_7.tmrca.50kb[4:ncol(rep_7.tmrca.50kb)], 1, mean)

rep_7.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_7.theta.50kb.bedgraph", header = T)
rep_7.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_7.rho.50kb.bedgraph", header = T)

rep_7.maps.50kb <- as.data.frame(cbind(rep_7.pi.50kb$sample_mean,
                                       rep_7.tmrca.50kb$sample_mean,
                                       rep_7.theta.50kb$sample_mean,
                                       rep_7.rho.50kb$sample_mean))

names(rep_7.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_7.maps.50kb$thetaC <- rep_7.maps.50kb$theta - mean(rep_7.maps.50kb$theta)
rep_7.maps.50kb$tmrcaC <- rep_7.maps.50kb$tmrca - mean(rep_7.maps.50kb$tmrca)
rep_7.maps.50kb$rhoC <- rep_7.maps.50kb$rho - mean(rep_7.maps.50kb$rho)

rep_7.maps.50kb$bin <- 1:nrow(rep_7.maps.50kb)

m.rep_7.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_7.maps.50kb)

rep_7.anova <- Anova(m.rep_7.50kb)
apiss <- rep_7.anova$"Sum Sq"
rep_7.anova$VarExp <- apiss / sum(apiss)




## rep_8

rep_8.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_8.diversity.50kb.bedgraph", header = T)
rep_8.pi.50kb$sample_mean <- apply(rep_8.pi.50kb[,(4:ncol(rep_8.pi.50kb))], 1, mean)

rep_8.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_8.TMRCA.50kb.bedgraph", header = T)
rep_8.tmrca.50kb$sample_mean <- apply(rep_8.tmrca.50kb[4:ncol(rep_8.tmrca.50kb)], 1, mean)

rep_8.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_8.theta.50kb.bedgraph", header = T)
rep_8.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_8.rho.50kb.bedgraph", header = T)

rep_8.maps.50kb <- as.data.frame(cbind(rep_8.pi.50kb$sample_mean,
                                       rep_8.tmrca.50kb$sample_mean,
                                       rep_8.theta.50kb$sample_mean,
                                       rep_8.rho.50kb$sample_mean))

names(rep_8.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_8.maps.50kb$thetaC <- rep_8.maps.50kb$theta - mean(rep_8.maps.50kb$theta)
rep_8.maps.50kb$tmrcaC <- rep_8.maps.50kb$tmrca - mean(rep_8.maps.50kb$tmrca)
rep_8.maps.50kb$rhoC <- rep_8.maps.50kb$rho - mean(rep_8.maps.50kb$rho)

rep_8.maps.50kb$bin <- 1:nrow(rep_8.maps.50kb)

m.rep_8.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_8.maps.50kb)

rep_8.anova <- Anova(m.rep_8.50kb)
apiss <- rep_8.anova$"Sum Sq"
rep_8.anova$VarExp <- apiss / sum(apiss)



## rep_9

rep_9.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_9.diversity.50kb.bedgraph", header = T)
rep_9.pi.50kb$sample_mean <- apply(rep_9.pi.50kb[,(4:ncol(rep_9.pi.50kb))], 1, mean)

rep_9.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_9.TMRCA.50kb.bedgraph", header = T)
rep_9.tmrca.50kb$sample_mean <- apply(rep_9.tmrca.50kb[4:ncol(rep_9.tmrca.50kb)], 1, mean)

rep_9.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_9.theta.50kb.bedgraph", header = T)
rep_9.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_9.rho.50kb.bedgraph", header = T)

rep_9.maps.50kb <- as.data.frame(cbind(rep_9.pi.50kb$sample_mean,
                                       rep_9.tmrca.50kb$sample_mean,
                                       rep_9.theta.50kb$sample_mean,
                                       rep_9.rho.50kb$sample_mean))

names(rep_9.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_9.maps.50kb$thetaC <- rep_9.maps.50kb$theta - mean(rep_9.maps.50kb$theta)
rep_9.maps.50kb$tmrcaC <- rep_9.maps.50kb$tmrca - mean(rep_9.maps.50kb$tmrca)
rep_9.maps.50kb$rhoC <- rep_9.maps.50kb$rho - mean(rep_9.maps.50kb$rho)

rep_9.maps.50kb$bin <- 1:nrow(rep_9.maps.50kb)

m.rep_9.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_9.maps.50kb)

rep_9.anova <- Anova(m.rep_9.50kb)
apiss <- rep_9.anova$"Sum Sq"
rep_9.anova$VarExp <- apiss / sum(apiss)


rep_4.anova


## rep_1

rep_1.pi.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_1.diversity.50kb.bedgraph", header = T)
rep_1.pi.50kb$sample_mean <- apply(rep_1.pi.50kb[,(4:ncol(rep_1.pi.50kb))], 1, mean)

rep_1.tmrca.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_1.TMRCA.50kb.bedgraph", header = T)
rep_1.tmrca.50kb$sample_mean <- apply(rep_1.tmrca.50kb[4:ncol(rep_1.tmrca.50kb)], 1, mean)

rep_1.theta.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_1.theta.50kb.bedgraph", header = T)
rep_1.rho.50kb <- read.table("flat_mu_inf_maps/dm_2L_bgs_rep_1.rho.50kb.bedgraph", header = T)

rep_1.maps.50kb <- as.data.frame(cbind(rep_1.pi.50kb$sample_mean,
                                       rep_1.tmrca.50kb$sample_mean,
                                       rep_1.theta.50kb$sample_mean,
                                       rep_1.rho.50kb$sample_mean))

names(rep_1.maps.50kb) <- c("diversity", "tmrca", "theta", "rho")

# centering
rep_1.maps.50kb$thetaC <- rep_1.maps.50kb$theta - mean(rep_1.maps.50kb$theta)
rep_1.maps.50kb$tmrcaC <- rep_1.maps.50kb$tmrca - mean(rep_1.maps.50kb$tmrca)
rep_1.maps.50kb$rhoC <- rep_1.maps.50kb$rho - mean(rep_1.maps.50kb$rho)

rep_1.maps.50kb$bin <- 1:nrow(rep_1.maps.50kb)

m.rep_1.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = rep_1.maps.50kb)

rep_1.anova <- Anova(m.rep_1.50kb)
apiss <- rep_1.anova$"Sum Sq"
rep_1.anova$VarExp <- apiss / sum(apiss)



rep_1.anova










# Plots
scale.5d <- function(x) sprintf("%.5f", x) # digits shown in y axis

molten.diversity <- melt(bgs.lands.50kb[c(1,8)], id.vars = "bin")
diversity.map.50kb <- ggplot(data = molten.diversity, aes(x = bin * 5e+4, y = value))
diversity.map.50kb <- diversity.map.50kb + geom_line(data = molten.diversity, colour = "#636363")
diversity.map.50kb <- diversity.map.50kb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
diversity.map.50kb <- diversity.map.50kb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
diversity.map.50kb <- diversity.map.50kb + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map.50kb <- diversity.map.50kb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 1e-3, yend = 1e-3, size = 2))
diversity.map.50kb <- diversity.map.50kb + theme(axis.title.x=element_blank(),
                                                 axis.text.x=element_blank(),
                                                 axis.ticks.x=element_blank(),
                                                 axis.title.y = element_text(size = 20), legend.position = "none") 

molten.rho <- melt(bgs.lands.50kb[c(3,8)], id.vars = "bin")
rho.map.50kb <- ggplot(data = molten.rho, aes(x = bin * 5e+4, y = value))
rho.map.50kb <- rho.map.50kb + geom_line(data = molten.rho, colour = "#636363")
rho.map.50kb <- rho.map.50kb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
rho.map.50kb <- rho.map.50kb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
rho.map.50kb <- rho.map.50kb + labs(title = NULL, x = NULL, y = expression(rho))
#rho.map.50kb <- rho.map.50kb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 2.2e-4, yend = 2.2e-4, size = 2))
rho.map.50kb <- rho.map.50kb + theme(axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank(),
                                     axis.title.y = element_text(size = 20), legend.position = "none") 

molten.theta <- melt(bgs.lands.50kb[c(2,8)], id.vars = "bin")
theta.map.50kb <- ggplot(data = molten.theta, aes(x = bin * 5e+4, y = value)) 
theta.map.50kb <- theta.map.50kb + geom_line(data = molten.theta, colour = "#636363")
theta.map.50kb <- theta.map.50kb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
theta.map.50kb <- theta.map.50kb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
theta.map.50kb <- theta.map.50kb + labs(title = NULL, x = NULL, y = expression(theta)) 
#theta.map.50kb <- theta.map.50kb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 1.45e-4, yend = 1.45e-4, size = 2))
theta.map.50kb <- theta.map.50kb + theme(axis.title.x=element_blank(),
                                         axis.text.x=element_blank(),
                                         axis.ticks.x=element_blank(), axis.title.y = element_text(size = 20), legend.position = "none") 

molten.tmrca <- melt(bgs.lands.50kb[c(4,8)], id.vars = "bin")
tmrca.map.50kb <- ggplot(data = molten.tmrca, aes(x = bin * 5e+4, y = value))
tmrca.map.50kb <- tmrca.map.50kb + geom_line(data = molten.tmrca, colour = "#636363")
tmrca.map.50kb <- tmrca.map.50kb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
tmrca.map.50kb <- tmrca.map.50kb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
tmrca.map.50kb <- tmrca.map.50kb + labs(title = NULL, x = "Position (bp)", y = expression(tau))
#tmrca.map.50kb <- tmrca.map.50kb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 4, yend = 4, size = 2))
tmrca.map.50kb <- tmrca.map.50kb + theme(axis.title.y = element_text(size = 20), legend.position = "none") 

maps.50kb <- plot_grid(diversity.map.50kb, rho.map.50kb, theta.map.50kb, tmrca.map.50kb, nrow = 4)
maps.50kb
save_plot("bgs.50kb.pdf", maps.50kb, base_height = 8, base_width = 16)


# 200kb
bgs.lands.200kb <- as.data.frame(cbind(bgs.pi.200kb$mean, bgs.theta.200kb$sample_mean, bgs.rho.200kb$sample_mean, bgs.TMRCA.200kb$mean))
names(bgs.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

# centering
bgs.lands.200kb$thetaC <- bgs.lands.200kb$theta - mean(bgs.lands.200kb$theta)
bgs.lands.200kb$tmrcaC <- bgs.lands.200kb$tmrca - mean(bgs.lands.200kb$tmrca)
bgs.lands.200kb$rhoC <- bgs.lands.200kb$rho - mean(bgs.lands.200kb$rho)

bgs.lands.200kb$bin <- 1:nrow(bgs.lands.200kb)

m.bgs.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = bgs.lands.200kb)
m.bgs.200kb.2 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC + rhoC:tmrcaC, data = bgs.lands.200kb)
m.bgs.200kb.3 <- lm(diversity ~ (thetaC + rhoC + tmrcaC) ^ 2, data = bgs.lands.200kb)

AIC(m.bgs.200kb, m.bgs.200kb.2, m.bgs.200kb.3)

bgs.lands.200kb.scaled <- bgs.lands.200kb * 100 # because otherwise precision issues arise in the ANOVA function
m.bgs.200kb.scaled <- lm(diversity ~ (thetaC + rhoC + tmrcaC) ^ 2, data = bgs.lands.200kb.scaled)

plot(resid(m.bgs.200kb.scaled)~fitted(m.bgs.200kb.scaled))
dwtest(m.bgs.200kb.scaled)
hmctest(m.bgs.200kb.scaled)
hist(resid(m.bgs.200kb.scaled))

summary(m.bgs.200kb.scaled)

anova.diversity <- Anova(m.bgs.200kb.scaled)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

#Anova Table (Type II tests)

#Response: diversity
#Sum Sq  Df   F value  Pr(>F)  VarExp
#thetaC        0.0001573   1  177.8794 0.00000 0.02311
#rhoC          0.0000001   1    0.1167 0.73329 0.00002
#tmrcaC        0.0065187   1 7370.0254 0.00000 0.95769
#thetaC:rhoC   0.0000001   1    0.0840 0.77248 0.00001
#thetaC:tmrcaC 0.0000311   1   35.1262 0.00000 0.00456
#rhoC:tmrcaC   0.0000013   1    1.4334 0.23376 0.00019
#Residuals     0.0000982 111                   0.01442

g.bgs.200kb.1 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = bgs.lands.200kb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

g.bgs.200kb.2 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = bgs.lands.200kb, weights = varPower(0, ~theta), cor = corAR1(0, ~bin), method = "ML")

g.bgs.200kb.3 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = bgs.lands.200kb, weights = varPower(0, ~theta), method = "ML")

g.bgs.200kb.4 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = bgs.lands.200kb, weights = varPower(0, ~tmrcaC), method = "ML")

AIC(g.bgs.200kb.1, g.bgs.200kb.2, g.bgs.200kb.3, g.bgs.200kb.4)

summary(g.bgs.200kb.2)

#Coefficients:
#  Value Std.Error   t-value p-value
#(Intercept)    0.000119 0.0000009 134.58700  0.0000
#thetaC         5.333327 0.5699165   9.35809  0.0000
#rhoC          -0.003620 0.0207195  -0.17474  0.8616
#tmrcaC         0.000147 0.0000017  84.22464  0.0000
#thetaC:tmrcaC  5.791087 0.8257040   7.01351  0.0000

vif(g.bgs.200kb.2)
#thetaC          rhoC        tmrcaC thetaC:tmrcaC 
#1.200503      1.010727      1.173875      1.287284

# Plots
molten.diversity <- melt(bgs.lands.200kb[c(1,8)], id.vars = "bin")
diversity.map.200kb <- ggplot(data = molten.diversity, aes(x = bin * 2e+5, y = value))
diversity.map.200kb <- diversity.map.200kb + geom_line(data = molten.diversity, colour = "#636363")
diversity.map.200kb <- diversity.map.200kb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
diversity.map.200kb <- diversity.map.200kb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
diversity.map.200kb <- diversity.map.200kb + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map.200kb <- diversity.map.200kb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 5e-4, yend = 5e-4, size = 2))
diversity.map.200kb <- diversity.map.200kb + theme(axis.title.x=element_blank(),
                                                 axis.text.x=element_blank(),
                                                 axis.ticks.x=element_blank(),
                                                 axis.title.y = element_text(size = 20), legend.position = "none") 

molten.rho <- melt(bgs.lands.200kb[c(3,8)], id.vars = "bin")
rho.map.200kb <- ggplot(data = molten.rho, aes(x = bin * 2e+5, y = value))
rho.map.200kb <- rho.map.200kb + geom_line(data = molten.rho, colour = "#636363")
rho.map.200kb <- rho.map.200kb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
rho.map.200kb <- rho.map.200kb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
rho.map.200kb <- rho.map.200kb + labs(title = NULL, x = NULL, y = expression(rho))
#rho.map.200kb <- rho.map.200kb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 2.2e-4, yend = 2.2e-4, size = 2))
rho.map.200kb <- rho.map.200kb + theme(axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank(),
                                     axis.title.y = element_text(size = 20), legend.position = "none") 

molten.theta <- melt(bgs.lands.200kb[c(2,8)], id.vars = "bin")
theta.map.200kb <- ggplot(data = molten.theta, aes(x = bin * 2e+5, y = value)) 
theta.map.200kb <- theta.map.200kb + geom_line(data = molten.theta, colour = "#636363")
theta.map.200kb <- theta.map.200kb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
theta.map.200kb <- theta.map.200kb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
theta.map.200kb <- theta.map.200kb + labs(title = NULL, x = NULL, y = expression(theta)) 
#theta.map.200kb <- theta.map.200kb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 1.45e-4, yend = 1.45e-4, size = 2))
theta.map.200kb <- theta.map.200kb + theme(axis.title.x=element_blank(),
                                         axis.text.x=element_blank(),
                                         axis.ticks.x=element_blank(), axis.title.y = element_text(size = 20), legend.position = "none") 

molten.tmrca <- melt(bgs.lands.200kb[c(4,8)], id.vars = "bin")
tmrca.map.200kb <- ggplot(data = molten.tmrca, aes(x = bin * 2e+5, y = value))
tmrca.map.200kb <- tmrca.map.200kb + geom_line(data = molten.tmrca, colour = "#636363")
tmrca.map.200kb <- tmrca.map.200kb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
tmrca.map.200kb <- tmrca.map.200kb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
tmrca.map.200kb <- tmrca.map.200kb + labs(title = NULL, x = "Position (bp)", y = expression(tau))
#tmrca.map.200kb <- tmrca.map.200kb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 4, yend = 4, size = 2))
tmrca.map.200kb <- tmrca.map.200kb + theme(axis.title.y = element_text(size = 20), legend.position = "none") 

maps.200kb <- plot_grid(diversity.map.200kb, rho.map.200kb, theta.map.200kb, tmrca.map.200kb, nrow = 4)
maps.200kb
save_plot("bgs.200kb.pdf", maps.200kb, base_height = 8, base_width = 16)



# 1Mb
bgs.lands.1Mb <- as.data.frame(cbind(bgs.pi.1Mb$mean, bgs.theta.1Mb$sample_mean, bgs.rho.1Mb$sample_mean, bgs.TMRCA.1Mb$mean))
names(bgs.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

# centering
bgs.lands.1Mb$thetaC <- bgs.lands.1Mb$theta - mean(bgs.lands.1Mb$theta)
bgs.lands.1Mb$tmrcaC <- bgs.lands.1Mb$tmrca - mean(bgs.lands.1Mb$tmrca)
bgs.lands.1Mb$rhoC <- bgs.lands.1Mb$rho - mean(bgs.lands.1Mb$rho)

bgs.lands.1Mb$bin <- 1:nrow(bgs.lands.1Mb)

m.bgs.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = bgs.lands.1Mb)
m.bgs.1Mb.2 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC + rhoC:tmrcaC, data = bgs.lands.1Mb)
m.bgs.1Mb.3 <- lm(diversity ~ (thetaC + rhoC + tmrcaC) ^ 2, data = bgs.lands.1Mb)

AIC(m.bgs.1Mb, m.bgs.1Mb.2, m.bgs.1Mb.3)

bgs.lands.1Mb.scaled <- bgs.lands.1Mb * 100 # because otherwise precision issues arise in the ANOVA function
m.bgs.1Mb.scaled <- lm(diversity ~ (thetaC + rhoC + tmrcaC) ^ 2, data = bgs.lands.1Mb.scaled)

plot(resid(m.bgs.1Mb.scaled)~fitted(m.bgs.1Mb.scaled))
dwtest(m.bgs.1Mb.scaled)
hmctest(m.bgs.1Mb.scaled)
hist(resid(m.bgs.1Mb.scaled))

summary(m.bgs.1Mb.scaled)

anova.diversity <- Anova(m.bgs.1Mb.scaled)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

#Anova Table (Type II tests)

#Response: diversity
#Sum Sq Df   F value  Pr(>F)  VarExp
#thetaC        1.1441e-05  1   56.8575 0.00000 0.03571
#rhoC          2.1900e-07  1    1.0904 0.31101 0.00068
#tmrcaC        3.0502e-04  1 1515.8585 0.00000 0.95198
#thetaC:rhoC   1.8100e-07  1    0.8975 0.35673 0.00056
#thetaC:tmrcaC 4.7000e-08  1    0.2319 0.63629 0.00015
#rhoC:tmrcaC   7.8000e-08  1    0.3853 0.54301 0.00024
#Residuals     3.4210e-06 17                   0.01068

g.bgs.1Mb.1 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                     data = bgs.lands.1Mb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

g.bgs.1Mb.2 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                     data = bgs.lands.1Mb, weights = varPower(0, ~theta), cor = corAR1(0, ~bin), method = "ML")

g.bgs.1Mb.3 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                     data = bgs.lands.1Mb, weights = varPower(0, ~theta), method = "ML")

g.bgs.1Mb.4 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                     data = bgs.lands.1Mb, weights = varPower(0, ~tmrcaC), method = "ML")

AIC(g.bgs.1Mb.1, g.bgs.1Mb.2, g.bgs.1Mb.3, g.bgs.1Mb.4)

summary(g.bgs.1Mb.4)

#Coefficients:
#  Value Std.Error   t-value p-value
#(Intercept)    0.000121  0.000001 132.33758  0.0000
#thetaC         7.198040  1.587574   4.53399  0.0002
#rhoC          -0.029352  0.027663  -1.06105  0.3020
#tmrcaC         0.000141  0.000004  37.98559  0.0000
#thetaC:tmrcaC  7.790249  4.792801   1.62541  0.1205

vif(g.bgs.1Mb.4)
#thetaC          rhoC        tmrcaC thetaC:tmrcaC 
#2.025202      1.179691      1.290478      1.882174 

# Plots
molten.diversity <- melt(bgs.lands.1Mb[c(1,8)], id.vars = "bin")
diversity.map.1Mb <- ggplot(data = molten.diversity, aes(x = bin * 1e+6, y = value))
diversity.map.1Mb <- diversity.map.1Mb + geom_line(data = molten.diversity, colour = "#636363")
diversity.map.1Mb <- diversity.map.1Mb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
diversity.map.1Mb <- diversity.map.1Mb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
diversity.map.1Mb <- diversity.map.1Mb + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map.1Mb <- diversity.map.1Mb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 2.25e-4, yend = 2.25e-4, size = 2))
diversity.map.1Mb <- diversity.map.1Mb + theme(axis.title.x=element_blank(),
                                                   axis.text.x=element_blank(),
                                                   axis.ticks.x=element_blank(),
                                                   axis.title.y = element_text(size = 20), legend.position = "none") 

molten.rho <- melt(bgs.lands.1Mb[c(3,8)], id.vars = "bin")
rho.map.1Mb <- ggplot(data = molten.rho, aes(x = bin * 1e+6, y = value))
rho.map.1Mb <- rho.map.1Mb + geom_line(data = molten.rho, colour = "#636363")
rho.map.1Mb <- rho.map.1Mb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
rho.map.1Mb <- rho.map.1Mb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
rho.map.1Mb <- rho.map.1Mb + labs(title = NULL, x = NULL, y = expression(rho))
#rho.map.1Mb <- rho.map.1Mb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 2.2e-4, yend = 2.2e-4, size = 2))
rho.map.1Mb <- rho.map.1Mb + theme(axis.title.x=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.ticks.x=element_blank(),
                                       axis.title.y = element_text(size = 20), legend.position = "none") 

molten.theta <- melt(bgs.lands.1Mb[c(2,8)], id.vars = "bin")
theta.map.1Mb <- ggplot(data = molten.theta, aes(x = bin * 1e+6, y = value)) 
theta.map.1Mb <- theta.map.1Mb + geom_line(data = molten.theta, colour = "#636363")
theta.map.1Mb <- theta.map.1Mb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
theta.map.1Mb <- theta.map.1Mb + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
theta.map.1Mb <- theta.map.1Mb + labs(title = NULL, x = NULL, y = expression(theta)) 
#theta.map.1Mb <- theta.map.1Mb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 1.45e-4, yend = 1.45e-4, size = 2))
theta.map.1Mb <- theta.map.1Mb + theme(axis.title.x=element_blank(),
                                           axis.text.x=element_blank(),
                                           axis.ticks.x=element_blank(), axis.title.y = element_text(size = 20), legend.position = "none") 

molten.tmrca <- melt(bgs.lands.1Mb[c(4,8)], id.vars = "bin")
tmrca.map.1Mb <- ggplot(data = molten.tmrca, aes(x = bin * 1e+6, y = value))
tmrca.map.1Mb <- tmrca.map.1Mb + geom_line(data = molten.tmrca, colour = "#636363")
tmrca.map.1Mb <- tmrca.map.1Mb + geom_smooth(method = "loess", se = F, colour = "#636363") + theme_bw()
tmrca.map.1Mb <- tmrca.map.1Mb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.5d)
tmrca.map.1Mb <- tmrca.map.1Mb + labs(title = NULL, x = "Position (bp)", y = expression(tau))
#tmrca.map.1Mb <- tmrca.map.1Mb + geom_segment(data = dm_2L_exome, aes(x = exon_start, xend = exon_end, y = 4, yend = 4, size = 2))
tmrca.map.1Mb <- tmrca.map.1Mb + theme(axis.title.y = element_text(size = 20), legend.position = "none") 

maps.1Mb <- plot_grid(diversity.map.1Mb, rho.map.1Mb, theta.map.1Mb, tmrca.map.1Mb, nrow = 4)
maps.1Mb
save_plot("bgs.1Mb.pdf", maps.1Mb, base_height = 8, base_width = 16)


