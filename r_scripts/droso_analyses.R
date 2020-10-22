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
library(RColorBrewer)

setwd("~")
setwd("Data/iSMC/theta_paper/real_data/")

# R_2 table for plotting at the end
r2.chr.tab <- as.data.frame(matrix(ncol = 5, nrow = 3))
colnames(r2.chr.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Bin_size(kb)")

theme.blank <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nreps <- 10

reps <- character(length = nreps)
for(i in 1:nreps) {
  reps[i] <- paste("rep_", i, sep = "")
}

##################c####################
#
# Drosophila-like neutral simulations of 2L vs Real Data 2L
#
########################################

# building linear models

# 50 kb
r2.inf.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 5))
row.names(r2.inf.50kb) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA")
colnames(r2.inf.50kb) <- reps

# sim landscapes
sim.rho.50k <- read.table("dm.sim.rho.50000.txt")
sim.theta.50k <- read.table("dm.sim.theta.50000.txt")
sim.lands.50k <- as.data.frame(cbind(sim.theta.50k$sim, sim.rho.50k$sim))
names(sim.lands.50k) <- c("theta", "rho")


# rep 1
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_1/rs.pair_1.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)
  
inf.lands.50k.rep1 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep1) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep1$thetaC <- inf.lands.50k.rep1$theta - mean(inf.lands.50k.rep1$theta)
inf.lands.50k.rep1$tmrcaC <- inf.lands.50k.rep1$tmrca - mean(inf.lands.50k.rep1$tmrca)
inf.lands.50k.rep1$rhoC <- inf.lands.50k.rep1$rho - mean(inf.lands.50k.rep1$rho)
  
inf.lands.50k.rep1$bin <- 1:nrow(inf.lands.50k.rep1)

cor.test(~theta+diversity, data = inf.lands.50k.rep1, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep1, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep1, method = "spearman")

m.diversity.rep_1 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.50k.rep1)

plot(resid(m.diversity.rep_1)~fitted(m.diversity.rep_1))
dwtest(m.diversity.rep_1)
hmctest(m.diversity.rep_1)
hist(resid(m.diversity.rep_1))
  
summary(m.diversity.rep_1) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.043e-02  2.236e-05 913.668   <2e-16 ***
# thetaC         1.096e+00  2.364e-03 463.564   <2e-16 ***
# rhoC          -4.461e-03  1.756e-02  -0.254      0.8    
# tmrcaC         2.026e-02  1.640e-04 123.483   <2e-16 ***
# thetaC:tmrcaC  1.100e+00  1.587e-02  69.290   <2e-16 ***

interact_plot(m.diversity.rep_1, pred = thetaC, modx= tmrcaC)

g.rep_1 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep1, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_1)
# (Intercept)   0.0204388 0.000031255 653.9355  0.0000
# thetaC        1.0980191 0.002985790 367.7482  0.0000
# tmrcaC        0.0199958 0.000176128 113.5300  0.0000
# rhoC          0.0035408 0.017040672   0.2078  0.8355
# thetaC:tmrcaC 1.0653107 0.017126226  62.2035  0.0000

anova.diversity <- Anova(m.diversity.rep_1)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)
  
r2.inf.50kb[1, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 1] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 1] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 1] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 1] <- anova.diversity$VarExp[4] * 100


# rep 2
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_2/rs.pair_2.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep2 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep2) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep2$thetaC <- inf.lands.50k.rep2$theta - mean(inf.lands.50k.rep2$theta)
inf.lands.50k.rep2$tmrcaC <- inf.lands.50k.rep2$tmrca - mean(inf.lands.50k.rep2$tmrca)
inf.lands.50k.rep2$rhoC <- inf.lands.50k.rep2$rho - mean(inf.lands.50k.rep2$rho)

inf.lands.50k.rep2$bin <- 1:nrow(inf.lands.50k.rep2)

cor.test(~theta+diversity, data = inf.lands.50k.rep2, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep2, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep2, method = "spearman")

m.diversity.rep_2 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.50k.rep2)

plot(resid(m.diversity.rep_2)~fitted(m.diversity.rep_2))
dwtest(m.diversity.rep_2)
hmctest(m.diversity.rep_2)
hist(resid(m.diversity.rep_2))

summary(m.diversity.rep_2) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.043e-02  2.069e-05 987.265   <2e-16 *** 
# thetaC         1.087e+00  2.167e-03 501.668   <2e-16 ***
# rhoC          -1.167e-02  1.648e-02  -0.708    0.479    
# tmrcaC         1.999e-02  1.568e-04 127.438   <2e-16 ***
# thetaC:tmrcaC  1.063e+00  1.365e-02  77.841   <2e-16 ***

interact_plot(m.diversity.rep_2, pred = thetaC, modx= tmrcaC)

g.rep_2 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep2, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_2)
# Estimate   Std.Error  t-value p-value
# (Intercept)    0.0204377 0.000030697 665.7833  0.0000
# thetaC         1.0884205 0.002800847 388.6041  0.0000
# tmrcaC         0.0196845 0.000163272 120.5625  0.0000
# rhoC          -0.0167357 0.015485360  -1.0807  0.2802
# thetaC:tmrcaC  1.0175254 0.014295984  71.1756  0.0000

anova.diversity <- Anova(m.diversity.rep_2)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 2] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 2] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 2] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 2] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 2] <- anova.diversity$VarExp[4] * 100
  

# rep_3
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_3/rs.pair_3.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep3 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep3) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep3$thetaC <- inf.lands.50k.rep3$theta - mean(inf.lands.50k.rep3$theta)
inf.lands.50k.rep3$tmrcaC <- inf.lands.50k.rep3$tmrca - mean(inf.lands.50k.rep3$tmrca)
inf.lands.50k.rep3$rhoC <- inf.lands.50k.rep3$rho - mean(inf.lands.50k.rep3$rho)

inf.lands.50k.rep3$bin <- 1:nrow(inf.lands.50k.rep3)

cor.test(~theta+diversity, data = inf.lands.50k.rep3, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep3, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep3, method = "spearman")

m.diversity.rep_3 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.50k.rep3)

plot(resid(m.diversity.rep_3)~fitted(m.diversity.rep_3))
dwtest(m.diversity.rep_3)
hmctest(m.diversity.rep_3)
hist(resid(m.diversity.rep_3))

summary(m.diversity.rep_3) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    0.0205593  0.0000198 1038.471   <2e-16 ***
# thetaC         1.0769710  0.0020324  529.902   <2e-16 ***
# rhoC          -0.0138881  0.0157297   -0.883    0.378    
# tmrcaC         0.0202409  0.0001408  143.768   <2e-16 ***
# thetaC:tmrcaC  1.0174596  0.0120993   84.092   <2e-16 ***

interact_plot(m.diversity.rep_3, pred = thetaC, modx = tmrcaC)
  
g.rep_3 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep3, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_3)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0205629 0.000026507 775.7537  0.0000
# thetaC         1.0781251 0.002506188 430.1853  0.0000
# tmrcaC         0.0201223 0.000150857 133.3870  0.0000
# rhoC          -0.0178210 0.015324110  -1.1629  0.2453
# thetaC:tmrcaC  0.9985038 0.013270574  75.2419  0.0000

anova.diversity <- Anova(m.diversity.rep_3)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 3] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 3] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 3] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 3] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 3] <- anova.diversity$VarExp[4] * 100


# rep_4
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_4/rs.pair_4.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep4 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep4) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep4$thetaC <- inf.lands.50k.rep4$theta - mean(inf.lands.50k.rep4$theta)
inf.lands.50k.rep4$tmrcaC <- inf.lands.50k.rep4$tmrca - mean(inf.lands.50k.rep4$tmrca)
inf.lands.50k.rep4$rhoC <- inf.lands.50k.rep4$rho - mean(inf.lands.50k.rep4$rho)

inf.lands.50k.rep4$bin <- 1:nrow(inf.lands.50k.rep4)

cor.test(~theta+diversity, data = inf.lands.50k.rep4, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep4, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep4, method = "spearman")

m.diversity.rep_4 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.50k.rep4)

plot(resid(m.diversity.rep_4)~fitted(m.diversity.rep_4))
dwtest(m.diversity.rep_4)
hmctest(m.diversity.rep_4)
hist(resid(m.diversity.rep_4))

summary(m.diversity.rep_4) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.061e-02  2.304e-05 894.744   <2e-16 ***
# thetaC         1.124e+00  2.418e-03 464.707   <2e-16 ***
# rhoC          -1.112e-02  1.845e-02  -0.603    0.547    
# tmrcaC         1.960e-02  1.576e-04 124.378   <2e-16 ***
# thetaC:tmrcaC  1.055e+00  1.499e-02  70.356   <2e-16 ***

interact_plot(m.diversity.rep_4, pred = thetaC, modx= tmrcaC)

g.rep_4 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep4, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_4)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0206153 0.000042516 484.8841  0.0000
# thetaC         1.1223093 0.003267770 343.4480  0.0000
# tmrcaC         0.0195383 0.000157223 124.2717  0.0000
# rhoC          -0.0073321 0.015064681  -0.4867  0.6266
# thetaC:tmrcaC  1.0070650 0.014339665  70.2293  0.0000

anova.diversity <- Anova(m.diversity.rep_4)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 4] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 4] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 4] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 4] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 4] <- anova.diversity$VarExp[4] * 100


# rep_5
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_5/rs.pair_5.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep5 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep5) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep5$thetaC <- inf.lands.50k.rep5$theta - mean(inf.lands.50k.rep5$theta)
inf.lands.50k.rep5$tmrcaC <- inf.lands.50k.rep5$tmrca - mean(inf.lands.50k.rep5$tmrca)
inf.lands.50k.rep5$rhoC <- inf.lands.50k.rep5$rho - mean(inf.lands.50k.rep5$rho)

inf.lands.50k.rep5$bin <- 1:nrow(inf.lands.50k.rep5)

cor.test(~theta+diversity, data = inf.lands.50k.rep5, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep5, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep5, method = "spearman")

m.diversity.rep_5 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.50k.rep5)

plot(resid(m.diversity.rep_5)~fitted(m.diversity.rep_5))
dwtest(m.diversity.rep_5)
hmctest(m.diversity.rep_5)
hist(resid(m.diversity.rep_5))

summary(m.diversity.rep_5) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    0.0204679  0.0000189 1082.787   <2e-16 ***
# thetaC         1.1078321  0.0019892  556.934   <2e-16 ***
# rhoC          -0.0069343  0.0156805   -0.442    0.658
# tmrcaC         0.0198380  0.0001375  144.318   <2e-16 ***
# thetaC:tmrcaC  1.0621447  0.0118730   89.459   <2e-16 ***

interact_plot(m.diversity.rep_5, pred = thetaC, modx= tmrcaC)

g.rep_5 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep5, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_5)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204733 0.000028264 724.3584  0.0000
# thetaC         1.1067703 0.002603042 425.1835  0.0000
# tmrcaC         0.0195333 0.000148186 131.8159  0.0000
# rhoC          -0.0117411 0.014941182  -0.7858  0.4323
# thetaC:tmrcaC  1.0119684 0.012760307  79.3060  0.0000

anova.diversity <- Anova(m.diversity.rep_5)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 5] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 5] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 5] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 5] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 5] <- anova.diversity$VarExp[4] * 100


# rep_6
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_6/rs.pair_6.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep6 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep6) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep6$thetaC <- inf.lands.50k.rep6$theta - mean(inf.lands.50k.rep6$theta)
inf.lands.50k.rep6$tmrcaC <- inf.lands.50k.rep6$tmrca - mean(inf.lands.50k.rep6$tmrca)
inf.lands.50k.rep6$rhoC <- inf.lands.50k.rep6$rho - mean(inf.lands.50k.rep6$rho)

inf.lands.50k.rep6$bin <- 1:nrow(inf.lands.50k.rep6)

cor.test(~theta+diversity, data = inf.lands.50k.rep6, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep6, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep6, method = "spearman")

m.diversity.rep_6 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.50k.rep6)

plot(resid(m.diversity.rep_6)~fitted(m.diversity.rep_6))
dwtest(m.diversity.rep_6)
hmctest(m.diversity.rep_6)
hist(resid(m.diversity.rep_6))

summary(m.diversity.rep_6) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.046e-02  2.183e-05 937.468   <2e-16 ***
# thetaC        1.086e+00  2.249e-03 482.709   <2e-16 ***
# rhoC          9.344e-03  1.682e-02   0.555    0.579    
# tmrcaC        2.016e-02  1.691e-04 119.259   <2e-16 ***
# thetaC:tmrcaC 1.050e+00  1.470e-02  71.423   <2e-16 ***

interact_plot(m.diversity.rep_6, pred = thetaC, modx= tmrcaC)

g.rep_6 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep6, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_6)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0204681 0.000030586 669.1903  0.0000
# thetaC        1.0858244 0.002832026 383.4091  0.0000
# tmrcaC        0.0199376 0.000176983 112.6524  0.0000
# rhoC          0.0002293 0.016283078   0.0141  0.9888
# thetaC:tmrcaC 1.0151180 0.015909552  63.8056  0.0000

anova.diversity <- Anova(m.diversity.rep_6)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 6] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 6] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 6] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 6] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 6] <- anova.diversity$VarExp[4] * 100



# rep_7
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_7/rs.pair_7.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep7 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep7) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep7$thetaC <- inf.lands.50k.rep7$theta - mean(inf.lands.50k.rep7$theta)
inf.lands.50k.rep7$tmrcaC <- inf.lands.50k.rep7$tmrca - mean(inf.lands.50k.rep7$tmrca)
inf.lands.50k.rep7$rhoC <- inf.lands.50k.rep7$rho - mean(inf.lands.50k.rep7$rho)

inf.lands.50k.rep7$bin <- 1:nrow(inf.lands.50k.rep7)

cor.test(~theta+diversity, data = inf.lands.50k.rep7, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep7, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep7, method = "spearman")

m.diversity.rep_7 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.50k.rep7)

plot(resid(m.diversity.rep_7)~fitted(m.diversity.rep_7))
dwtest(m.diversity.rep_7)
hmctest(m.diversity.rep_7)
hist(resid(m.diversity.rep_7))

summary(m.diversity.rep_7) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    2.046e-02  1.906e-05 1073.378  < 2e-16 ***
# thetaC         1.093e+00  1.994e-03  548.323  < 2e-16 ***
# rhoC          -3.984e-02  1.517e-02   -2.627  0.00885 ** 
# tmrcaC         2.021e-02  1.462e-04  138.265  < 2e-16 ***
# thetaC:tmrcaC  1.056e+00  1.344e-02   78.551  < 2e-16 ***

interact_plot(m.diversity.rep_7, pred = thetaC, modx= tmrcaC)

g.rep_7 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep7, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_7)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204701 0.000032313 633.4967  0.0000
# thetaC         1.0964625 0.002720215 403.0793  0.0000
# tmrcaC         0.0197764 0.000149360 132.4076  0.0000
# rhoC          -0.0310546 0.013470566  -2.3054  0.0215
# thetaC:tmrcaC  0.9959236 0.014228586  69.9946  0.0000

anova.diversity <- Anova(m.diversity.rep_7)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 7] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 7] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 7] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 7] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 7] <- anova.diversity$VarExp[4] * 100



# rep_8
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_8/rs.pair_8.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep8 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep8) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep8$thetaC <- inf.lands.50k.rep8$theta - mean(inf.lands.50k.rep8$theta)
inf.lands.50k.rep8$tmrcaC <- inf.lands.50k.rep8$tmrca - mean(inf.lands.50k.rep8$tmrca)
inf.lands.50k.rep8$rhoC <- inf.lands.50k.rep8$rho - mean(inf.lands.50k.rep8$rho)

inf.lands.50k.rep8$bin <- 1:nrow(inf.lands.50k.rep8)


cor.test(~theta+diversity, data = inf.lands.50k.rep8, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep8, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep8, method = "spearman")

m.diversity.rep_8 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.50k.rep8)

plot(resid(m.diversity.rep_8)~fitted(m.diversity.rep_8))
dwtest(m.diversity.rep_8)
hmctest(m.diversity.rep_8)
hist(resid(m.diversity.rep_8))

summary(m.diversity.rep_8) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.070e-02  1.938e-05 1067.98   <2e-16 ***
# thetaC         1.107e+00  2.023e-03  547.27   <2e-16 ***
# rhoC          -6.776e-03  1.539e-02   -0.44     0.66    
# tmrcaC         2.013e-02  1.393e-04  144.50   <2e-16 ***
# thetaC:tmrcaC  1.071e+00  1.222e-02   87.61   <2e-16 ***

interact_plot(m.diversity.rep_8, pred = thetaC, modx= tmrcaC)

g.rep_8 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep8, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_8)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0207031 0.000029299 706.6176  0.0000
# thetaC        1.1082548 0.002641688 419.5253  0.0000
# tmrcaC        0.0198182 0.000148295 133.6401  0.0000
# rhoC          0.0047060 0.014548476   0.3235  0.7465
# thetaC:tmrcaC 1.0223058 0.013460362  75.9494  0.0000

anova.diversity <- Anova(m.diversity.rep_8)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 8] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 8] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 8] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 8] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 8] <- anova.diversity$VarExp[4] * 100



# rep_9
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_9/rs.pair_9.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep9 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep9) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep9$thetaC <- inf.lands.50k.rep9$theta - mean(inf.lands.50k.rep9$theta)
inf.lands.50k.rep9$tmrcaC <- inf.lands.50k.rep9$tmrca - mean(inf.lands.50k.rep9$tmrca)
inf.lands.50k.rep9$rhoC <- inf.lands.50k.rep9$rho - mean(inf.lands.50k.rep9$rho)

inf.lands.50k.rep9$bin <- 1:nrow(inf.lands.50k.rep9)

cor.test(~theta+diversity, data = inf.lands.50k.rep9, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep9, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep9, method = "spearman")

m.diversity.rep_9 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.50k.rep9)

plot(resid(m.diversity.rep_9)~fitted(m.diversity.rep_9))
dwtest(m.diversity.rep_9)
hmctest(m.diversity.rep_9)
hist(resid(m.diversity.rep_9))

summary(m.diversity.rep_9) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.037e-02  2.069e-05 984.748   <2e-16 ***
# thetaC        1.086e+00  2.147e-03 505.918   <2e-16 ***
# rhoC          1.171e-02  1.718e-02   0.682    0.496    
# tmrcaC        1.990e-02  1.617e-04 123.090   <2e-16 ***
 #thetaC:tmrcaC 1.045e+00  1.408e-02  74.239   <2e-16 ***

interact_plot(m.diversity.rep_9, pred = thetaC, modx= tmrcaC)

g.rep_9 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep9, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_9)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0203835 0.000032571 625.8123  0.0000
# thetaC        1.0864837 0.002863160 379.4701  0.0000
# tmrcaC        0.0194098 0.000164202 118.2072  0.0000
# rhoC          0.0014542 0.015946130   0.0912  0.9274
# thetaC:tmrcaC 0.9860233 0.014806396  66.5944  0.0000

anova.diversity <- Anova(m.diversity.rep_9)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 9] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 9] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 9] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 9] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 9] <- anova.diversity$VarExp[4] * 100

# rep_10
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_10/rs.pair_10.", sep = "")
pi.50k <- read.table(paste(p, "diversity.50kb.bedgraph", sep = ""), header = T)
pi.50k$avg <- apply(pi.50k[4:ncol(pi.50k)], 1, mean)
tmrca.50k <- read.table(paste(p, "TMRCA.50kb.bedgraph", sep = ""), header = T)
tmrca.50k$avg <- apply(tmrca.50k[4:ncol(tmrca.50k)], 1, mean)
rho.50k <- read.table(paste(p, "rho.50kb.bedgraph", sep = ""), header = T)
theta.50k <- read.table(paste(p, "theta.50kb.bedgraph", sep = ""), header = T)

inf.lands.50k.rep10 <- as.data.frame(cbind(pi.50k$avg, theta.50k$sample_mean, rho.50k$sample_mean, tmrca.50k$avg))
names(inf.lands.50k.rep10) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.50k.rep10$thetaC <- inf.lands.50k.rep10$theta - mean(inf.lands.50k.rep10$theta)
inf.lands.50k.rep10$tmrcaC <- inf.lands.50k.rep10$tmrca - mean(inf.lands.50k.rep10$tmrca)
inf.lands.50k.rep10$rhoC <- inf.lands.50k.rep10$rho - mean(inf.lands.50k.rep10$rho)

inf.lands.50k.rep10$bin <- 1:nrow(inf.lands.50k.rep10)

cor.test(~theta+diversity, data = inf.lands.50k.rep10, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.50k.rep10, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.50k.rep10, method = "spearman")

m.diversity.rep_10 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.50k.rep10)

plot(resid(m.diversity.rep_10)~fitted(m.diversity.rep_10))
dwtest(m.diversity.rep_10)
hmctest(m.diversity.rep_10)
hist(resid(m.diversity.rep_10))

summary(m.diversity.rep_10) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.048e-02  2.313e-05 885.243   <2e-16 ***
# thetaC         1.074e+00  2.407e-03 446.391   <2e-16 ***
# rhoC          -2.139e-02  1.771e-02  -1.208    0.227    
# tmrcaC         2.070e-02  1.767e-04 117.152   <2e-16 ***
# thetaC:tmrcaC  1.096e+00  1.562e-02  70.181   <2e-16 ***

interact_plot(m.diversity.rep_10, pred = thetaC, modx= tmrcaC)

g.rep_10 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.50k.rep10, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_10)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204923 0.000036424 562.6011  0.0000
# thetaC         1.0766673 0.003238378 332.4711  0.0000
# tmrcaC         0.0201003 0.000188420 106.6779  0.0000
# rhoC          -0.0294997 0.016448791  -1.7934  0.0734
# thetaC:tmrcaC  1.0263521 0.016424787  62.4880  0.0000

anova.diversity <- Anova(m.diversity.rep_10)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.50kb[1, 10] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.50kb[2, 10] <- anova.diversity$VarExp[1] * 100
r2.inf.50kb[3, 10] <- anova.diversity$VarExp[2] * 100
r2.inf.50kb[4, 10] <- anova.diversity$VarExp[3] * 100
r2.inf.50kb[5, 10] <- anova.diversity$VarExp[4] * 100

r2.inf.50kb$average <- rowMeans(r2.inf.50kb)
r2.inf.50kb <- transform(r2.inf.50kb, sd=apply(r2.inf.50kb, 1, sd, na.rm = TRUE))

# plots
rho.plot <- as.data.frame(cbind(inf.lands.50k.rep1$bin,
                                sim.rho.50k$sim, 
                                inf.lands.50k.rep1$rho,
                                inf.lands.50k.rep2$rho,
                                inf.lands.50k.rep3$rho,
                                inf.lands.50k.rep4$rho,
                                inf.lands.50k.rep5$rho,
                                inf.lands.50k.rep6$rho,
                                inf.lands.50k.rep7$rho,
                                inf.lands.50k.rep8$rho,
                                inf.lands.50k.rep9$rho,
                                inf.lands.50k.rep10$rho))

rho.plot[,3:ncol(rho.plot)] <- 1.53 * rho.plot[,3:ncol(rho.plot)]

names(rho.plot) <- c("bin", "sim", reps)
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map.50kb <- ggplot(data = molten.rho, aes(x = bin * 50, y = value, colour = variable)) + theme.blank
rho.map.50kb <- rho.map.50kb + geom_line(data = molten.rho, aes(size = variable)) + scale_size_manual(values = c(1.2, rep(0.4, 10)))
rho.map.50kb <- rho.map.50kb + scale_color_manual(values = c("black", brewer.pal(n = 9, name = "Blues"), "cyan"))
rho.map.50kb <- rho.map.50kb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map.50kb <- rho.map.50kb + labs(title = NULL, x = NULL, y = expression(rho))
rho.map.50kb <- rho.map.50kb + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

theta.plot <- as.data.frame(cbind(inf.lands.50k.rep1$bin,
                                sim.theta.50k$sim, 
                                inf.lands.50k.rep1$theta,
                                inf.lands.50k.rep2$theta,
                                inf.lands.50k.rep3$theta,
                                inf.lands.50k.rep4$theta,
                                inf.lands.50k.rep5$theta,
                                inf.lands.50k.rep6$theta,
                                inf.lands.50k.rep7$theta,
                                inf.lands.50k.rep8$theta,
                                inf.lands.50k.rep9$theta,
                                inf.lands.50k.rep10$theta))

names(theta.plot) <- c("bin", "sim", reps)
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map.50kb <- ggplot(data = molten.theta, aes(x = bin * 50, y = value, colour = variable)) + theme.blank
theta.map.50kb <- theta.map.50kb + geom_line(data = molten.theta, aes(size = variable)) + scale_size_manual(values = c(1.2, rep(0.4, 10)))
theta.map.50kb <- theta.map.50kb + scale_color_manual(values = c("black", brewer.pal(n = 9, name = "Reds"), "orange"))
theta.map.50kb <- theta.map.50kb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map.50kb <- theta.map.50kb + labs(title = NULL, x = NULL, y = expression(theta))
theta.map.50kb <- theta.map.50kb + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))




# 200 kb
r2.inf.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 5))
row.names(r2.inf.200kb) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA")
colnames(r2.inf.200kb) <- reps

# sim landscapes
sim.rho.200k <- read.table("dm.sim.rho.2e+05.txt")
sim.theta.200k <- read.table("dm.sim.theta.2e+05.txt")
sim.lands.200k <- as.data.frame(cbind(sim.theta.200k$sim, sim.rho.200k$sim))
names(sim.lands.200k) <- c("theta", "rho")


# rep 1
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_1/rs.pair_1.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep1 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep1) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep1$thetaC <- inf.lands.200k.rep1$theta - mean(inf.lands.200k.rep1$theta)
inf.lands.200k.rep1$tmrcaC <- inf.lands.200k.rep1$tmrca - mean(inf.lands.200k.rep1$tmrca)
inf.lands.200k.rep1$rhoC <- inf.lands.200k.rep1$rho - mean(inf.lands.200k.rep1$rho)

inf.lands.200k.rep1$bin <- 1:nrow(inf.lands.200k.rep1)

cor.test(~theta+diversity, data = inf.lands.200k.rep1, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep1, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep1, method = "spearman")

m.diversity.rep_1 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.200k.rep1)

plot(resid(m.diversity.rep_1)~fitted(m.diversity.rep_1))
dwtest(m.diversity.rep_1)
hmctest(m.diversity.rep_1)
hist(resid(m.diversity.rep_1))

summary(m.diversity.rep_1) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.043e-02  2.236e-05 913.668   <2e-16 ***
# thetaC         1.096e+00  2.364e-03 463.564   <2e-16 ***
# rhoC          -4.461e-03  1.756e-02  -0.254      0.8    
# tmrcaC         2.026e-02  1.640e-04 123.483   <2e-16 ***
# thetaC:tmrcaC  1.100e+00  1.587e-02  69.290   <2e-16 ***

interact_plot(m.diversity.rep_1, pred = thetaC, modx= tmrcaC)

g.rep_1 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep1, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_1)
# (Intercept)   0.0204388 0.000031255 653.9355  0.0000
# thetaC        1.0980191 0.002985790 367.7482  0.0000
# tmrcaC        0.0199958 0.000176128 113.5300  0.0000
# rhoC          0.0035408 0.017040672   0.2078  0.8355
# thetaC:tmrcaC 1.0653107 0.017126226  62.2035  0.0000

anova.diversity <- Anova(m.diversity.rep_1)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 1] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 1] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 1] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 1] <- anova.diversity$VarExp[4] * 100


# rep 2
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_2/rs.pair_2.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep2 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep2) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep2$thetaC <- inf.lands.200k.rep2$theta - mean(inf.lands.200k.rep2$theta)
inf.lands.200k.rep2$tmrcaC <- inf.lands.200k.rep2$tmrca - mean(inf.lands.200k.rep2$tmrca)
inf.lands.200k.rep2$rhoC <- inf.lands.200k.rep2$rho - mean(inf.lands.200k.rep2$rho)

inf.lands.200k.rep2$bin <- 1:nrow(inf.lands.200k.rep2)

cor.test(~theta+diversity, data = inf.lands.200k.rep2, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep2, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep2, method = "spearman")

m.diversity.rep_2 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.200k.rep2)

plot(resid(m.diversity.rep_2)~fitted(m.diversity.rep_2))
dwtest(m.diversity.rep_2)
hmctest(m.diversity.rep_2)
hist(resid(m.diversity.rep_2))

summary(m.diversity.rep_2) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.043e-02  2.069e-05 987.265   <2e-16 *** 
# thetaC         1.087e+00  2.167e-03 501.668   <2e-16 ***
# rhoC          -1.167e-02  1.648e-02  -0.708    0.479    
# tmrcaC         1.999e-02  1.568e-04 127.438   <2e-16 ***
# thetaC:tmrcaC  1.063e+00  1.365e-02  77.841   <2e-16 ***

interact_plot(m.diversity.rep_2, pred = thetaC, modx= tmrcaC)

g.rep_2 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep2, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_2)
# Estimate   Std.Error  t-value p-value
# (Intercept)    0.0204377 0.000030697 665.7833  0.0000
# thetaC         1.0884205 0.002800847 388.6041  0.0000
# tmrcaC         0.0196845 0.000163272 120.5625  0.0000
# rhoC          -0.0167357 0.015485360  -1.0807  0.2802
# thetaC:tmrcaC  1.0175254 0.014295984  71.1756  0.0000

anova.diversity <- Anova(m.diversity.rep_2)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 2] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 2] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 2] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 2] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 2] <- anova.diversity$VarExp[4] * 100


# rep_3
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_3/rs.pair_3.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep3 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep3) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep3$thetaC <- inf.lands.200k.rep3$theta - mean(inf.lands.200k.rep3$theta)
inf.lands.200k.rep3$tmrcaC <- inf.lands.200k.rep3$tmrca - mean(inf.lands.200k.rep3$tmrca)
inf.lands.200k.rep3$rhoC <- inf.lands.200k.rep3$rho - mean(inf.lands.200k.rep3$rho)

inf.lands.200k.rep3$bin <- 1:nrow(inf.lands.200k.rep3)

cor.test(~theta+diversity, data = inf.lands.200k.rep3, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep3, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep3, method = "spearman")

m.diversity.rep_3 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.200k.rep3)

plot(resid(m.diversity.rep_3)~fitted(m.diversity.rep_3))
dwtest(m.diversity.rep_3)
hmctest(m.diversity.rep_3)
hist(resid(m.diversity.rep_3))

summary(m.diversity.rep_3) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    0.0205593  0.0000198 1038.471   <2e-16 ***
# thetaC         1.0769710  0.0020324  529.902   <2e-16 ***
# rhoC          -0.0138881  0.0157297   -0.883    0.378    
# tmrcaC         0.0202409  0.0001408  143.768   <2e-16 ***
# thetaC:tmrcaC  1.0174596  0.0120993   84.092   <2e-16 ***

interact_plot(m.diversity.rep_3, pred = thetaC, modx = tmrcaC)

g.rep_3 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep3, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_3)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0205629 0.000026507 775.7537  0.0000
# thetaC         1.0781251 0.002506188 430.1853  0.0000
# tmrcaC         0.0201223 0.000150857 133.3870  0.0000
# rhoC          -0.0178210 0.015324110  -1.1629  0.2453
# thetaC:tmrcaC  0.9985038 0.013270574  75.2419  0.0000

anova.diversity <- Anova(m.diversity.rep_3)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 3] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 3] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 3] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 3] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 3] <- anova.diversity$VarExp[4] * 100


# rep_4
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_4/rs.pair_4.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep4 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep4) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep4$thetaC <- inf.lands.200k.rep4$theta - mean(inf.lands.200k.rep4$theta)
inf.lands.200k.rep4$tmrcaC <- inf.lands.200k.rep4$tmrca - mean(inf.lands.200k.rep4$tmrca)
inf.lands.200k.rep4$rhoC <- inf.lands.200k.rep4$rho - mean(inf.lands.200k.rep4$rho)

inf.lands.200k.rep4$bin <- 1:nrow(inf.lands.200k.rep4)

cor.test(~theta+diversity, data = inf.lands.200k.rep4, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep4, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep4, method = "spearman")

m.diversity.rep_4 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.200k.rep4)

plot(resid(m.diversity.rep_4)~fitted(m.diversity.rep_4))
dwtest(m.diversity.rep_4)
hmctest(m.diversity.rep_4)
hist(resid(m.diversity.rep_4))

summary(m.diversity.rep_4) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.061e-02  2.304e-05 894.744   <2e-16 ***
# thetaC         1.124e+00  2.418e-03 464.707   <2e-16 ***
# rhoC          -1.112e-02  1.845e-02  -0.603    0.547    
# tmrcaC         1.960e-02  1.576e-04 124.378   <2e-16 ***
# thetaC:tmrcaC  1.055e+00  1.499e-02  70.356   <2e-16 ***

interact_plot(m.diversity.rep_4, pred = thetaC, modx= tmrcaC)

g.rep_4 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep4, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_4)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0206153 0.000042516 484.8841  0.0000
# thetaC         1.1223093 0.003267770 343.4480  0.0000
# tmrcaC         0.0195383 0.000157223 124.2717  0.0000
# rhoC          -0.0073321 0.015064681  -0.4867  0.6266
# thetaC:tmrcaC  1.0070650 0.014339665  70.2293  0.0000

anova.diversity <- Anova(m.diversity.rep_4)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 4] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 4] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 4] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 4] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 4] <- anova.diversity$VarExp[4] * 100


# rep_5
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_5/rs.pair_5.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep5 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep5) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep5$thetaC <- inf.lands.200k.rep5$theta - mean(inf.lands.200k.rep5$theta)
inf.lands.200k.rep5$tmrcaC <- inf.lands.200k.rep5$tmrca - mean(inf.lands.200k.rep5$tmrca)
inf.lands.200k.rep5$rhoC <- inf.lands.200k.rep5$rho - mean(inf.lands.200k.rep5$rho)

inf.lands.200k.rep5$bin <- 1:nrow(inf.lands.200k.rep5)

cor.test(~theta+diversity, data = inf.lands.200k.rep5, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep5, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep5, method = "spearman")

m.diversity.rep_5 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.200k.rep5)

plot(resid(m.diversity.rep_5)~fitted(m.diversity.rep_5))
dwtest(m.diversity.rep_5)
hmctest(m.diversity.rep_5)
hist(resid(m.diversity.rep_5))

summary(m.diversity.rep_5) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    0.0204679  0.0000189 1082.787   <2e-16 ***
# thetaC         1.1078321  0.0019892  556.934   <2e-16 ***
# rhoC          -0.0069343  0.0156805   -0.442    0.658
# tmrcaC         0.0198380  0.0001375  144.318   <2e-16 ***
# thetaC:tmrcaC  1.0621447  0.0118730   89.459   <2e-16 ***

interact_plot(m.diversity.rep_5, pred = thetaC, modx= tmrcaC)

g.rep_5 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep5, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_5)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204733 0.000028264 724.3584  0.0000
# thetaC         1.1067703 0.002603042 425.1835  0.0000
# tmrcaC         0.0195333 0.000148186 131.8159  0.0000
# rhoC          -0.0117411 0.014941182  -0.7858  0.4323
# thetaC:tmrcaC  1.0119684 0.012760307  79.3060  0.0000

anova.diversity <- Anova(m.diversity.rep_5)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 5] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 5] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 5] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 5] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 5] <- anova.diversity$VarExp[4] * 100


# rep_6
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_6/rs.pair_6.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep6 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep6) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep6$thetaC <- inf.lands.200k.rep6$theta - mean(inf.lands.200k.rep6$theta)
inf.lands.200k.rep6$tmrcaC <- inf.lands.200k.rep6$tmrca - mean(inf.lands.200k.rep6$tmrca)
inf.lands.200k.rep6$rhoC <- inf.lands.200k.rep6$rho - mean(inf.lands.200k.rep6$rho)

inf.lands.200k.rep6$bin <- 1:nrow(inf.lands.200k.rep6)

cor.test(~theta+diversity, data = inf.lands.200k.rep6, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep6, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep6, method = "spearman")

m.diversity.rep_6 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.200k.rep6)

plot(resid(m.diversity.rep_6)~fitted(m.diversity.rep_6))
dwtest(m.diversity.rep_6)
hmctest(m.diversity.rep_6)
hist(resid(m.diversity.rep_6))

summary(m.diversity.rep_6) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.046e-02  2.183e-05 937.468   <2e-16 ***
# thetaC        1.086e+00  2.249e-03 482.709   <2e-16 ***
# rhoC          9.344e-03  1.682e-02   0.555    0.579    
# tmrcaC        2.016e-02  1.691e-04 119.259   <2e-16 ***
# thetaC:tmrcaC 1.050e+00  1.470e-02  71.423   <2e-16 ***

interact_plot(m.diversity.rep_6, pred = thetaC, modx= tmrcaC)

g.rep_6 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep6, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_6)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0204681 0.000030586 669.1903  0.0000
# thetaC        1.0858244 0.002832026 383.4091  0.0000
# tmrcaC        0.0199376 0.000176983 112.6524  0.0000
# rhoC          0.0002293 0.016283078   0.0141  0.9888
# thetaC:tmrcaC 1.0151180 0.015909552  63.8056  0.0000

anova.diversity <- Anova(m.diversity.rep_6)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 6] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 6] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 6] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 6] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 6] <- anova.diversity$VarExp[4] * 100



# rep_7
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_7/rs.pair_7.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep7 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep7) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep7$thetaC <- inf.lands.200k.rep7$theta - mean(inf.lands.200k.rep7$theta)
inf.lands.200k.rep7$tmrcaC <- inf.lands.200k.rep7$tmrca - mean(inf.lands.200k.rep7$tmrca)
inf.lands.200k.rep7$rhoC <- inf.lands.200k.rep7$rho - mean(inf.lands.200k.rep7$rho)

inf.lands.200k.rep7$bin <- 1:nrow(inf.lands.200k.rep7)

cor.test(~theta+diversity, data = inf.lands.200k.rep7, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep7, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep7, method = "spearman")

m.diversity.rep_7 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.200k.rep7)

plot(resid(m.diversity.rep_7)~fitted(m.diversity.rep_7))
dwtest(m.diversity.rep_7)
hmctest(m.diversity.rep_7)
hist(resid(m.diversity.rep_7))

summary(m.diversity.rep_7) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    2.046e-02  1.906e-05 1073.378  < 2e-16 ***
# thetaC         1.093e+00  1.994e-03  548.323  < 2e-16 ***
# rhoC          -3.984e-02  1.517e-02   -2.627  0.00885 ** 
# tmrcaC         2.021e-02  1.462e-04  138.265  < 2e-16 ***
# thetaC:tmrcaC  1.056e+00  1.344e-02   78.551  < 2e-16 ***

interact_plot(m.diversity.rep_7, pred = thetaC, modx= tmrcaC)

g.rep_7 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep7, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_7)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204701 0.000032313 633.4967  0.0000
# thetaC         1.0964625 0.002720215 403.0793  0.0000
# tmrcaC         0.0197764 0.000149360 132.4076  0.0000
# rhoC          -0.0310546 0.013470566  -2.3054  0.0215
# thetaC:tmrcaC  0.9959236 0.014228586  69.9946  0.0000

anova.diversity <- Anova(m.diversity.rep_7)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 7] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 7] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 7] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 7] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 7] <- anova.diversity$VarExp[4] * 100



# rep_8
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_8/rs.pair_8.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep8 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep8) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep8$thetaC <- inf.lands.200k.rep8$theta - mean(inf.lands.200k.rep8$theta)
inf.lands.200k.rep8$tmrcaC <- inf.lands.200k.rep8$tmrca - mean(inf.lands.200k.rep8$tmrca)
inf.lands.200k.rep8$rhoC <- inf.lands.200k.rep8$rho - mean(inf.lands.200k.rep8$rho)

inf.lands.200k.rep8$bin <- 1:nrow(inf.lands.200k.rep8)


cor.test(~theta+diversity, data = inf.lands.200k.rep8, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep8, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep8, method = "spearman")

m.diversity.rep_8 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.200k.rep8)

plot(resid(m.diversity.rep_8)~fitted(m.diversity.rep_8))
dwtest(m.diversity.rep_8)
hmctest(m.diversity.rep_8)
hist(resid(m.diversity.rep_8))

summary(m.diversity.rep_8) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.070e-02  1.938e-05 1067.98   <2e-16 ***
# thetaC         1.107e+00  2.023e-03  547.27   <2e-16 ***
# rhoC          -6.776e-03  1.539e-02   -0.44     0.66    
# tmrcaC         2.013e-02  1.393e-04  144.50   <2e-16 ***
# thetaC:tmrcaC  1.071e+00  1.222e-02   87.61   <2e-16 ***

interact_plot(m.diversity.rep_8, pred = thetaC, modx= tmrcaC)

g.rep_8 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep8, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_8)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0207031 0.000029299 706.6176  0.0000
# thetaC        1.1082548 0.002641688 419.5253  0.0000
# tmrcaC        0.0198182 0.000148295 133.6401  0.0000
# rhoC          0.0047060 0.014548476   0.3235  0.7465
# thetaC:tmrcaC 1.0223058 0.013460362  75.9494  0.0000

anova.diversity <- Anova(m.diversity.rep_8)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 8] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 8] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 8] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 8] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 8] <- anova.diversity$VarExp[4] * 100



# rep_9
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_9/rs.pair_9.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep9 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep9) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep9$thetaC <- inf.lands.200k.rep9$theta - mean(inf.lands.200k.rep9$theta)
inf.lands.200k.rep9$tmrcaC <- inf.lands.200k.rep9$tmrca - mean(inf.lands.200k.rep9$tmrca)
inf.lands.200k.rep9$rhoC <- inf.lands.200k.rep9$rho - mean(inf.lands.200k.rep9$rho)

inf.lands.200k.rep9$bin <- 1:nrow(inf.lands.200k.rep9)

cor.test(~theta+diversity, data = inf.lands.200k.rep9, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep9, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep9, method = "spearman")

m.diversity.rep_9 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.200k.rep9)

plot(resid(m.diversity.rep_9)~fitted(m.diversity.rep_9))
dwtest(m.diversity.rep_9)
hmctest(m.diversity.rep_9)
hist(resid(m.diversity.rep_9))

summary(m.diversity.rep_9) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.037e-02  2.069e-05 984.748   <2e-16 ***
# thetaC        1.086e+00  2.147e-03 505.918   <2e-16 ***
# rhoC          1.171e-02  1.718e-02   0.682    0.496    
# tmrcaC        1.990e-02  1.617e-04 123.090   <2e-16 ***
#thetaC:tmrcaC 1.045e+00  1.408e-02  74.239   <2e-16 ***

interact_plot(m.diversity.rep_9, pred = thetaC, modx= tmrcaC)

g.rep_9 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.200k.rep9, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_9)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0203835 0.000032571 625.8123  0.0000
# thetaC        1.0864837 0.002863160 379.4701  0.0000
# tmrcaC        0.0194098 0.000164202 118.2072  0.0000
# rhoC          0.0014542 0.015946130   0.0912  0.9274
# thetaC:tmrcaC 0.9860233 0.014806396  66.5944  0.0000

anova.diversity <- Anova(m.diversity.rep_9)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 9] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 9] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 9] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 9] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 9] <- anova.diversity$VarExp[4] * 100

# rep_10
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_10/rs.pair_10.", sep = "")
pi.200k <- read.table(paste(p, "diversity.200kb.bedgraph", sep = ""), header = T)
pi.200k$avg <- apply(pi.200k[4:ncol(pi.200k)], 1, mean)
tmrca.200k <- read.table(paste(p, "TMRCA.200kb.bedgraph", sep = ""), header = T)
tmrca.200k$avg <- apply(tmrca.200k[4:ncol(tmrca.200k)], 1, mean)
rho.200k <- read.table(paste(p, "rho.200kb.bedgraph", sep = ""), header = T)
theta.200k <- read.table(paste(p, "theta.200kb.bedgraph", sep = ""), header = T)

inf.lands.200k.rep10 <- as.data.frame(cbind(pi.200k$avg, theta.200k$sample_mean, rho.200k$sample_mean, tmrca.200k$avg))
names(inf.lands.200k.rep10) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.200k.rep10$thetaC <- inf.lands.200k.rep10$theta - mean(inf.lands.200k.rep10$theta)
inf.lands.200k.rep10$tmrcaC <- inf.lands.200k.rep10$tmrca - mean(inf.lands.200k.rep10$tmrca)
inf.lands.200k.rep10$rhoC <- inf.lands.200k.rep10$rho - mean(inf.lands.200k.rep10$rho)

inf.lands.200k.rep10$bin <- 1:nrow(inf.lands.200k.rep10)

cor.test(~theta+diversity, data = inf.lands.200k.rep10, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.200k.rep10, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.200k.rep10, method = "spearman")

m.diversity.rep_10 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.200k.rep10)

plot(resid(m.diversity.rep_10)~fitted(m.diversity.rep_10))
dwtest(m.diversity.rep_10)
hmctest(m.diversity.rep_10)
hist(resid(m.diversity.rep_10))

summary(m.diversity.rep_10) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.048e-02  2.313e-05 885.243   <2e-16 ***
# thetaC         1.074e+00  2.407e-03 446.391   <2e-16 ***
# rhoC          -2.139e-02  1.771e-02  -1.208    0.227    
# tmrcaC         2.070e-02  1.767e-04 117.152   <2e-16 ***
# thetaC:tmrcaC  1.096e+00  1.562e-02  70.181   <2e-16 ***

interact_plot(m.diversity.rep_10, pred = thetaC, modx= tmrcaC)

g.rep_10 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                data = inf.lands.200k.rep10, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_10)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204923 0.000036424 562.6011  0.0000
# thetaC         1.0766673 0.003238378 332.4711  0.0000
# tmrcaC         0.0201003 0.000188420 106.6779  0.0000
# rhoC          -0.0294997 0.016448791  -1.7934  0.0734
# thetaC:tmrcaC  1.0263521 0.016424787  62.4880  0.0000

anova.diversity <- Anova(m.diversity.rep_10)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.200kb[1, 10] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.200kb[2, 10] <- anova.diversity$VarExp[1] * 100
r2.inf.200kb[3, 10] <- anova.diversity$VarExp[2] * 100
r2.inf.200kb[4, 10] <- anova.diversity$VarExp[3] * 100
r2.inf.200kb[5, 10] <- anova.diversity$VarExp[4] * 100

r2.inf.200kb$average <- rowMeans(r2.inf.200kb)
r2.inf.200kb <- transform(r2.inf.200kb, sd=apply(r2.inf.200kb, 1, sd, na.rm = TRUE))

# plots
rho.plot <- as.data.frame(cbind(inf.lands.200k.rep1$bin,
                                sim.rho.200k$sim, 
                                inf.lands.200k.rep1$rho,
                                inf.lands.200k.rep2$rho,
                                inf.lands.200k.rep3$rho,
                                inf.lands.200k.rep4$rho,
                                inf.lands.200k.rep5$rho,
                                inf.lands.200k.rep6$rho,
                                inf.lands.200k.rep7$rho,
                                inf.lands.200k.rep8$rho,
                                inf.lands.200k.rep9$rho,
                                inf.lands.200k.rep10$rho))

rho.plot[,3:ncol(rho.plot)] <- 1.53 * rho.plot[,3:ncol(rho.plot)]

names(rho.plot) <- c("bin", "sim", reps)
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map.200kb <- ggplot(data = molten.rho, aes(x = bin * 200, y = value, colour = variable)) + theme.blank
rho.map.200kb <- rho.map.200kb + geom_line(data = molten.rho, aes(size = variable)) + scale_size_manual(values = c(1.2, rep(0.4, 10)))
rho.map.200kb <- rho.map.200kb + scale_color_manual(values = c("black", brewer.pal(n = 9, name = "Blues"), "cyan"))
rho.map.200kb <- rho.map.200kb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map.200kb <- rho.map.200kb + labs(title = NULL, x = NULL, y = expression(rho))
rho.map.200kb <- rho.map.200kb + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

theta.plot <- as.data.frame(cbind(inf.lands.200k.rep1$bin,
                                  sim.theta.200k$sim, 
                                  inf.lands.200k.rep1$theta,
                                  inf.lands.200k.rep2$theta,
                                  inf.lands.200k.rep3$theta,
                                  inf.lands.200k.rep4$theta,
                                  inf.lands.200k.rep5$theta,
                                  inf.lands.200k.rep6$theta,
                                  inf.lands.200k.rep7$theta,
                                  inf.lands.200k.rep8$theta,
                                  inf.lands.200k.rep9$theta,
                                  inf.lands.200k.rep10$theta))

names(theta.plot) <- c("bin", "sim", reps)
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map.200kb <- ggplot(data = molten.theta, aes(x = bin * 200, y = value, colour = variable)) + theme.blank
theta.map.200kb <- theta.map.200kb + geom_line(data = molten.theta, aes(size = variable)) + scale_size_manual(values = c(1.2, rep(0.4, 10)))
theta.map.200kb <- theta.map.200kb + scale_color_manual(values = c("black", brewer.pal(n = 9, name = "Reds"), "orange"))
theta.map.200kb <- theta.map.200kb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map.200kb <- theta.map.200kb + labs(title = NULL, x = NULL, y = expression(theta))
theta.map.200kb <- theta.map.200kb + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))




# 1 Mb
r2.inf.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 5))
row.names(r2.inf.1Mb) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA")
colnames(r2.inf.1Mb) <- reps

# sim landscapes
sim.rho.1M <- read.table("dm.sim.rho.1e+06.txt")
sim.theta.1M <- read.table("dm.sim.theta.1e+06.txt")
sim.lands.1M <- as.data.frame(cbind(sim.theta.1M$sim, sim.rho.1M$sim))
names(sim.lands.1M) <- c("theta", "rho")


# rep 1
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_1/rs.pair_1.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep1 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep1) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep1$thetaC <- inf.lands.1M.rep1$theta - mean(inf.lands.1M.rep1$theta)
inf.lands.1M.rep1$tmrcaC <- inf.lands.1M.rep1$tmrca - mean(inf.lands.1M.rep1$tmrca)
inf.lands.1M.rep1$rhoC <- inf.lands.1M.rep1$rho - mean(inf.lands.1M.rep1$rho)

inf.lands.1M.rep1$bin <- 1:nrow(inf.lands.1M.rep1)

cor.test(~theta+diversity, data = inf.lands.1M.rep1, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep1, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep1, method = "spearman")

m.diversity.rep_1 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.1M.rep1)

plot(resid(m.diversity.rep_1)~fitted(m.diversity.rep_1))
dwtest(m.diversity.rep_1)
hmctest(m.diversity.rep_1)
hist(resid(m.diversity.rep_1))

summary(m.diversity.rep_1) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.043e-02  2.236e-05 913.668   <2e-16 ***
# thetaC         1.096e+00  2.364e-03 463.564   <2e-16 ***
# rhoC          -4.461e-03  1.756e-02  -0.254      0.8    
# tmrcaC         2.026e-02  1.640e-04 123.483   <2e-16 ***
# thetaC:tmrcaC  1.100e+00  1.587e-02  69.290   <2e-16 ***

interact_plot(m.diversity.rep_1, pred = thetaC, modx= tmrcaC)

g.rep_1 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep1, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_1)
# (Intercept)   0.0204388 0.000031255 653.9355  0.0000
# thetaC        1.0980191 0.002985790 367.7482  0.0000
# tmrcaC        0.0199958 0.000176128 113.5300  0.0000
# rhoC          0.0035408 0.017040672   0.2078  0.8355
# thetaC:tmrcaC 1.0653107 0.017126226  62.2035  0.0000

anova.diversity <- Anova(m.diversity.rep_1)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 1] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 1] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 1] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 1] <- anova.diversity$VarExp[4] * 100


# rep 2
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_2/rs.pair_2.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep2 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep2) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep2$thetaC <- inf.lands.1M.rep2$theta - mean(inf.lands.1M.rep2$theta)
inf.lands.1M.rep2$tmrcaC <- inf.lands.1M.rep2$tmrca - mean(inf.lands.1M.rep2$tmrca)
inf.lands.1M.rep2$rhoC <- inf.lands.1M.rep2$rho - mean(inf.lands.1M.rep2$rho)

inf.lands.1M.rep2$bin <- 1:nrow(inf.lands.1M.rep2)

cor.test(~theta+diversity, data = inf.lands.1M.rep2, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep2, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep2, method = "spearman")

m.diversity.rep_2 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.1M.rep2)

plot(resid(m.diversity.rep_2)~fitted(m.diversity.rep_2))
dwtest(m.diversity.rep_2)
hmctest(m.diversity.rep_2)
hist(resid(m.diversity.rep_2))

summary(m.diversity.rep_2) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.043e-02  2.069e-05 987.265   <2e-16 *** 
# thetaC         1.087e+00  2.167e-03 501.668   <2e-16 ***
# rhoC          -1.167e-02  1.648e-02  -0.708    0.479    
# tmrcaC         1.999e-02  1.568e-04 127.438   <2e-16 ***
# thetaC:tmrcaC  1.063e+00  1.365e-02  77.841   <2e-16 ***

interact_plot(m.diversity.rep_2, pred = thetaC, modx= tmrcaC)

g.rep_2 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep2, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_2)
# Estimate   Std.Error  t-value p-value
# (Intercept)    0.0204377 0.000030697 665.7833  0.0000
# thetaC         1.0884205 0.002800847 388.6041  0.0000
# tmrcaC         0.0196845 0.000163272 120.5625  0.0000
# rhoC          -0.0167357 0.015485360  -1.0807  0.2802
# thetaC:tmrcaC  1.0175254 0.014295984  71.1756  0.0000

anova.diversity <- Anova(m.diversity.rep_2)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 2] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 2] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 2] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 2] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 2] <- anova.diversity$VarExp[4] * 100


# rep_3
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_3/rs.pair_3.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep3 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep3) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep3$thetaC <- inf.lands.1M.rep3$theta - mean(inf.lands.1M.rep3$theta)
inf.lands.1M.rep3$tmrcaC <- inf.lands.1M.rep3$tmrca - mean(inf.lands.1M.rep3$tmrca)
inf.lands.1M.rep3$rhoC <- inf.lands.1M.rep3$rho - mean(inf.lands.1M.rep3$rho)

inf.lands.1M.rep3$bin <- 1:nrow(inf.lands.1M.rep3)

cor.test(~theta+diversity, data = inf.lands.1M.rep3, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep3, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep3, method = "spearman")

m.diversity.rep_3 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.1M.rep3)

plot(resid(m.diversity.rep_3)~fitted(m.diversity.rep_3))
dwtest(m.diversity.rep_3)
hmctest(m.diversity.rep_3)
hist(resid(m.diversity.rep_3))

summary(m.diversity.rep_3) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    0.0205593  0.0000198 1038.471   <2e-16 ***
# thetaC         1.0769710  0.0020324  529.902   <2e-16 ***
# rhoC          -0.0138881  0.0157297   -0.883    0.378    
# tmrcaC         0.0202409  0.0001408  143.768   <2e-16 ***
# thetaC:tmrcaC  1.0174596  0.0120993   84.092   <2e-16 ***

interact_plot(m.diversity.rep_3, pred = thetaC, modx = tmrcaC)

g.rep_3 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep3, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_3)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0205629 0.000026507 775.7537  0.0000
# thetaC         1.0781251 0.002506188 430.1853  0.0000
# tmrcaC         0.0201223 0.000150857 133.3870  0.0000
# rhoC          -0.0178210 0.015324110  -1.1629  0.2453
# thetaC:tmrcaC  0.9985038 0.013270574  75.2419  0.0000

anova.diversity <- Anova(m.diversity.rep_3)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 3] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 3] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 3] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 3] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 3] <- anova.diversity$VarExp[4] * 100


# rep_4
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_4/rs.pair_4.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep4 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep4) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep4$thetaC <- inf.lands.1M.rep4$theta - mean(inf.lands.1M.rep4$theta)
inf.lands.1M.rep4$tmrcaC <- inf.lands.1M.rep4$tmrca - mean(inf.lands.1M.rep4$tmrca)
inf.lands.1M.rep4$rhoC <- inf.lands.1M.rep4$rho - mean(inf.lands.1M.rep4$rho)

inf.lands.1M.rep4$bin <- 1:nrow(inf.lands.1M.rep4)

cor.test(~theta+diversity, data = inf.lands.1M.rep4, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep4, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep4, method = "spearman")

m.diversity.rep_4 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.1M.rep4)

plot(resid(m.diversity.rep_4)~fitted(m.diversity.rep_4))
dwtest(m.diversity.rep_4)
hmctest(m.diversity.rep_4)
hist(resid(m.diversity.rep_4))

summary(m.diversity.rep_4) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.061e-02  2.304e-05 894.744   <2e-16 ***
# thetaC         1.124e+00  2.418e-03 464.707   <2e-16 ***
# rhoC          -1.112e-02  1.845e-02  -0.603    0.547    
# tmrcaC         1.960e-02  1.576e-04 124.378   <2e-16 ***
# thetaC:tmrcaC  1.055e+00  1.499e-02  70.356   <2e-16 ***

interact_plot(m.diversity.rep_4, pred = thetaC, modx= tmrcaC)

g.rep_4 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep4, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_4)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0206153 0.000042516 484.8841  0.0000
# thetaC         1.1223093 0.003267770 343.4480  0.0000
# tmrcaC         0.0195383 0.000157223 124.2717  0.0000
# rhoC          -0.0073321 0.015064681  -0.4867  0.6266
# thetaC:tmrcaC  1.0070650 0.014339665  70.2293  0.0000

anova.diversity <- Anova(m.diversity.rep_4)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 4] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 4] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 4] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 4] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 4] <- anova.diversity$VarExp[4] * 100


# rep_5
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_5/rs.pair_5.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep5 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep5) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep5$thetaC <- inf.lands.1M.rep5$theta - mean(inf.lands.1M.rep5$theta)
inf.lands.1M.rep5$tmrcaC <- inf.lands.1M.rep5$tmrca - mean(inf.lands.1M.rep5$tmrca)
inf.lands.1M.rep5$rhoC <- inf.lands.1M.rep5$rho - mean(inf.lands.1M.rep5$rho)

inf.lands.1M.rep5$bin <- 1:nrow(inf.lands.1M.rep5)

cor.test(~theta+diversity, data = inf.lands.1M.rep5, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep5, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep5, method = "spearman")

m.diversity.rep_5 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC, data = inf.lands.1M.rep5)

plot(resid(m.diversity.rep_5)~fitted(m.diversity.rep_5))
dwtest(m.diversity.rep_5)
hmctest(m.diversity.rep_5)
hist(resid(m.diversity.rep_5))

summary(m.diversity.rep_5) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    0.0204679  0.0000189 1082.787   <2e-16 ***
# thetaC         1.1078321  0.0019892  556.934   <2e-16 ***
# rhoC          -0.0069343  0.0156805   -0.442    0.658
# tmrcaC         0.0198380  0.0001375  144.318   <2e-16 ***
# thetaC:tmrcaC  1.0621447  0.0118730   89.459   <2e-16 ***

interact_plot(m.diversity.rep_5, pred = thetaC, modx= tmrcaC)

g.rep_5 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep5, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_5)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204733 0.000028264 724.3584  0.0000
# thetaC         1.1067703 0.002603042 425.1835  0.0000
# tmrcaC         0.0195333 0.000148186 131.8159  0.0000
# rhoC          -0.0117411 0.014941182  -0.7858  0.4323
# thetaC:tmrcaC  1.0119684 0.012760307  79.3060  0.0000

anova.diversity <- Anova(m.diversity.rep_5)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 5] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 5] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 5] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 5] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 5] <- anova.diversity$VarExp[4] * 100


# rep_6
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_6/rs.pair_6.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep6 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep6) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep6$thetaC <- inf.lands.1M.rep6$theta - mean(inf.lands.1M.rep6$theta)
inf.lands.1M.rep6$tmrcaC <- inf.lands.1M.rep6$tmrca - mean(inf.lands.1M.rep6$tmrca)
inf.lands.1M.rep6$rhoC <- inf.lands.1M.rep6$rho - mean(inf.lands.1M.rep6$rho)

inf.lands.1M.rep6$bin <- 1:nrow(inf.lands.1M.rep6)

cor.test(~theta+diversity, data = inf.lands.1M.rep6, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep6, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep6, method = "spearman")

m.diversity.rep_6 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.1M.rep6)

plot(resid(m.diversity.rep_6)~fitted(m.diversity.rep_6))
dwtest(m.diversity.rep_6)
hmctest(m.diversity.rep_6)
hist(resid(m.diversity.rep_6))

summary(m.diversity.rep_6) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.046e-02  2.183e-05 937.468   <2e-16 ***
# thetaC        1.086e+00  2.249e-03 482.709   <2e-16 ***
# rhoC          9.344e-03  1.682e-02   0.555    0.579    
# tmrcaC        2.016e-02  1.691e-04 119.259   <2e-16 ***
# thetaC:tmrcaC 1.050e+00  1.470e-02  71.423   <2e-16 ***

interact_plot(m.diversity.rep_6, pred = thetaC, modx= tmrcaC)

g.rep_6 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep6, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_6)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0204681 0.000030586 669.1903  0.0000
# thetaC        1.0858244 0.002832026 383.4091  0.0000
# tmrcaC        0.0199376 0.000176983 112.6524  0.0000
# rhoC          0.0002293 0.016283078   0.0141  0.9888
# thetaC:tmrcaC 1.0151180 0.015909552  63.8056  0.0000

anova.diversity <- Anova(m.diversity.rep_6)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 6] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 6] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 6] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 6] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 6] <- anova.diversity$VarExp[4] * 100



# rep_7
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_7/rs.pair_7.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep7 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep7) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep7$thetaC <- inf.lands.1M.rep7$theta - mean(inf.lands.1M.rep7$theta)
inf.lands.1M.rep7$tmrcaC <- inf.lands.1M.rep7$tmrca - mean(inf.lands.1M.rep7$tmrca)
inf.lands.1M.rep7$rhoC <- inf.lands.1M.rep7$rho - mean(inf.lands.1M.rep7$rho)

inf.lands.1M.rep7$bin <- 1:nrow(inf.lands.1M.rep7)

cor.test(~theta+diversity, data = inf.lands.1M.rep7, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep7, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep7, method = "spearman")

m.diversity.rep_7 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.1M.rep7)

plot(resid(m.diversity.rep_7)~fitted(m.diversity.rep_7))
dwtest(m.diversity.rep_7)
hmctest(m.diversity.rep_7)
hist(resid(m.diversity.rep_7))

summary(m.diversity.rep_7) 
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    2.046e-02  1.906e-05 1073.378  < 2e-16 ***
# thetaC         1.093e+00  1.994e-03  548.323  < 2e-16 ***
# rhoC          -3.984e-02  1.517e-02   -2.627  0.00885 ** 
# tmrcaC         2.021e-02  1.462e-04  138.265  < 2e-16 ***
# thetaC:tmrcaC  1.056e+00  1.344e-02   78.551  < 2e-16 ***

interact_plot(m.diversity.rep_7, pred = thetaC, modx= tmrcaC)

g.rep_7 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep7, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_7)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204701 0.000032313 633.4967  0.0000
# thetaC         1.0964625 0.002720215 403.0793  0.0000
# tmrcaC         0.0197764 0.000149360 132.4076  0.0000
# rhoC          -0.0310546 0.013470566  -2.3054  0.0215
# thetaC:tmrcaC  0.9959236 0.014228586  69.9946  0.0000

anova.diversity <- Anova(m.diversity.rep_7)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 7] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 7] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 7] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 7] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 7] <- anova.diversity$VarExp[4] * 100



# rep_8
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_8/rs.pair_8.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep8 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep8) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep8$thetaC <- inf.lands.1M.rep8$theta - mean(inf.lands.1M.rep8$theta)
inf.lands.1M.rep8$tmrcaC <- inf.lands.1M.rep8$tmrca - mean(inf.lands.1M.rep8$tmrca)
inf.lands.1M.rep8$rhoC <- inf.lands.1M.rep8$rho - mean(inf.lands.1M.rep8$rho)

inf.lands.1M.rep8$bin <- 1:nrow(inf.lands.1M.rep8)


cor.test(~theta+diversity, data = inf.lands.1M.rep8, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep8, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep8, method = "spearman")

m.diversity.rep_8 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.1M.rep8)

plot(resid(m.diversity.rep_8)~fitted(m.diversity.rep_8))
dwtest(m.diversity.rep_8)
hmctest(m.diversity.rep_8)
hist(resid(m.diversity.rep_8))

summary(m.diversity.rep_8) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.070e-02  1.938e-05 1067.98   <2e-16 ***
# thetaC         1.107e+00  2.023e-03  547.27   <2e-16 ***
# rhoC          -6.776e-03  1.539e-02   -0.44     0.66    
# tmrcaC         2.013e-02  1.393e-04  144.50   <2e-16 ***
# thetaC:tmrcaC  1.071e+00  1.222e-02   87.61   <2e-16 ***

interact_plot(m.diversity.rep_8, pred = thetaC, modx= tmrcaC)

g.rep_8 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep8, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_8)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0207031 0.000029299 706.6176  0.0000
# thetaC        1.1082548 0.002641688 419.5253  0.0000
# tmrcaC        0.0198182 0.000148295 133.6401  0.0000
# rhoC          0.0047060 0.014548476   0.3235  0.7465
# thetaC:tmrcaC 1.0223058 0.013460362  75.9494  0.0000

anova.diversity <- Anova(m.diversity.rep_8)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 8] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 8] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 8] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 8] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 8] <- anova.diversity$VarExp[4] * 100



# rep_9
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_9/rs.pair_9.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep9 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep9) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep9$thetaC <- inf.lands.1M.rep9$theta - mean(inf.lands.1M.rep9$theta)
inf.lands.1M.rep9$tmrcaC <- inf.lands.1M.rep9$tmrca - mean(inf.lands.1M.rep9$tmrca)
inf.lands.1M.rep9$rhoC <- inf.lands.1M.rep9$rho - mean(inf.lands.1M.rep9$rho)

inf.lands.1M.rep9$bin <- 1:nrow(inf.lands.1M.rep9)

cor.test(~theta+diversity, data = inf.lands.1M.rep9, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep9, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep9, method = "spearman")

m.diversity.rep_9 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.1M.rep9)

plot(resid(m.diversity.rep_9)~fitted(m.diversity.rep_9))
dwtest(m.diversity.rep_9)
hmctest(m.diversity.rep_9)
hist(resid(m.diversity.rep_9))

summary(m.diversity.rep_9) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.037e-02  2.069e-05 984.748   <2e-16 ***
# thetaC        1.086e+00  2.147e-03 505.918   <2e-16 ***
# rhoC          1.171e-02  1.718e-02   0.682    0.496    
# tmrcaC        1.990e-02  1.617e-04 123.090   <2e-16 ***
#thetaC:tmrcaC 1.045e+00  1.408e-02  74.239   <2e-16 ***

interact_plot(m.diversity.rep_9, pred = thetaC, modx= tmrcaC)

g.rep_9 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
               data = inf.lands.1M.rep9, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_9)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0203835 0.000032571 625.8123  0.0000
# thetaC        1.0864837 0.002863160 379.4701  0.0000
# tmrcaC        0.0194098 0.000164202 118.2072  0.0000
# rhoC          0.0014542 0.015946130   0.0912  0.9274
# thetaC:tmrcaC 0.9860233 0.014806396  66.5944  0.0000

anova.diversity <- Anova(m.diversity.rep_9)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 9] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 9] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 9] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 9] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 9] <- anova.diversity$VarExp[4] * 100

# rep_10
p <- paste("dm_chr_maps/2L/dm_coal_sims/rep_10/rs.pair_10.", sep = "")
pi.1M <- read.table(paste(p, "diversity.1Mb.bedgraph", sep = ""), header = T)
pi.1M$avg <- apply(pi.1M[4:ncol(pi.1M)], 1, mean)
tmrca.1M <- read.table(paste(p, "TMRCA.1Mb.bedgraph", sep = ""), header = T)
tmrca.1M$avg <- apply(tmrca.1M[4:ncol(tmrca.1M)], 1, mean)
rho.1M <- read.table(paste(p, "rho.1Mb.bedgraph", sep = ""), header = T)
theta.1M <- read.table(paste(p, "theta.1Mb.bedgraph", sep = ""), header = T)

inf.lands.1M.rep10 <- as.data.frame(cbind(pi.1M$avg, theta.1M$sample_mean, rho.1M$sample_mean, tmrca.1M$avg))
names(inf.lands.1M.rep10) <- c("diversity", "theta", "rho", "tmrca")

# centering
inf.lands.1M.rep10$thetaC <- inf.lands.1M.rep10$theta - mean(inf.lands.1M.rep10$theta)
inf.lands.1M.rep10$tmrcaC <- inf.lands.1M.rep10$tmrca - mean(inf.lands.1M.rep10$tmrca)
inf.lands.1M.rep10$rhoC <- inf.lands.1M.rep10$rho - mean(inf.lands.1M.rep10$rho)

inf.lands.1M.rep10$bin <- 1:nrow(inf.lands.1M.rep10)

cor.test(~theta+diversity, data = inf.lands.1M.rep10, method = "spearman")
cor.test(~tmrca+diversity, data = inf.lands.1M.rep10, method = "spearman")
cor.test(~rho+diversity, data = inf.lands.1M.rep10, method = "spearman")

m.diversity.rep_10 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC*tmrcaC , data = inf.lands.1M.rep10)

plot(resid(m.diversity.rep_10)~fitted(m.diversity.rep_10))
dwtest(m.diversity.rep_10)
hmctest(m.diversity.rep_10)
hist(resid(m.diversity.rep_10))

summary(m.diversity.rep_10) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.048e-02  2.313e-05 885.243   <2e-16 ***
# thetaC         1.074e+00  2.407e-03 446.391   <2e-16 ***
# rhoC          -2.139e-02  1.771e-02  -1.208    0.227    
# tmrcaC         2.070e-02  1.767e-04 117.152   <2e-16 ***
# thetaC:tmrcaC  1.096e+00  1.562e-02  70.181   <2e-16 ***

interact_plot(m.diversity.rep_10, pred = thetaC, modx= tmrcaC)

g.rep_10 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                data = inf.lands.1M.rep10, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_10)
# Value   Std.Error  t-value p-value
# (Intercept)    0.0204923 0.000036424 562.6011  0.0000
# thetaC         1.0766673 0.003238378 332.4711  0.0000
# tmrcaC         0.0201003 0.000188420 106.6779  0.0000
# rhoC          -0.0294997 0.016448791  -1.7934  0.0734
# thetaC:tmrcaC  1.0263521 0.016424787  62.4880  0.0000

anova.diversity <- Anova(m.diversity.rep_10)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.inf.1Mb[1, 10] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.inf.1Mb[2, 10] <- anova.diversity$VarExp[1] * 100
r2.inf.1Mb[3, 10] <- anova.diversity$VarExp[2] * 100
r2.inf.1Mb[4, 10] <- anova.diversity$VarExp[3] * 100
r2.inf.1Mb[5, 10] <- anova.diversity$VarExp[4] * 100

r2.inf.1Mb$average <- rowMeans(r2.inf.1Mb)
r2.inf.1Mb <- transform(r2.inf.1Mb, sd=apply(r2.inf.1Mb, 1, sd, na.rm = TRUE))

# plots
rho.plot <- as.data.frame(cbind(inf.lands.1M.rep1$bin,
                                sim.rho.1M$sim, 
                                 inf.lands.1M.rep1$rho,
                                 inf.lands.1M.rep2$rho,
                                 inf.lands.1M.rep3$rho,
                                 inf.lands.1M.rep4$rho,
                                 inf.lands.1M.rep5$rho,
                                 inf.lands.1M.rep6$rho,
                                 inf.lands.1M.rep7$rho,
                                 inf.lands.1M.rep8$rho,
                                 inf.lands.1M.rep9$rho,
                                 inf.lands.1M.rep10$rho))

rho.plot[,3:ncol(rho.plot)] <- 1.53 * rho.plot[,3:ncol(rho.plot)]

names(rho.plot) <- c("bin", "sim", reps)
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map.1Mb <- ggplot(data = molten.rho, aes(x = bin * 1000, y = value, colour = variable)) + theme.blank
rho.map.1Mb <- rho.map.1Mb + geom_line(data = molten.rho, aes(size = variable)) + scale_size_manual(values = c(1.2, rep(0.4, 10)))
rho.map.1Mb <- rho.map.1Mb + scale_color_manual(values = c("black", brewer.pal(n = 9, name = "Blues"), "cyan"))
rho.map.1Mb <- rho.map.1Mb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map.1Mb <- rho.map.1Mb + labs(title = NULL, x = NULL, y = expression(rho))
rho.map.1Mb <- rho.map.1Mb + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

theta.plot <- as.data.frame(cbind(inf.lands.1M.rep1$bin,
                                  sim.theta.1M$sim, 
                                  inf.lands.1M.rep1$theta,
                                  inf.lands.1M.rep2$theta,
                                  inf.lands.1M.rep3$theta,
                                  inf.lands.1M.rep4$theta,
                                  inf.lands.1M.rep5$theta,
                                  inf.lands.1M.rep6$theta,
                                  inf.lands.1M.rep7$theta,
                                  inf.lands.1M.rep8$theta,
                                  inf.lands.1M.rep9$theta,
                                  inf.lands.1M.rep10$theta))

names(theta.plot) <- c("bin", "sim", reps)
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map.1Mb <- ggplot(data = molten.theta, aes(x = bin * 1000, y = value, colour = variable)) + theme.blank
theta.map.1Mb <- theta.map.1Mb + geom_line(data = molten.theta, aes(size = variable)) + scale_size_manual(values = c(1.2, rep(0.4, 10)))
theta.map.1Mb <- theta.map.1Mb + scale_color_manual(values = c("black", brewer.pal(n = 9, name = "Reds"), "orange"))
theta.map.1Mb <- theta.map.1Mb + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map.1Mb <- theta.map.1Mb + labs(title = NULL, x = NULL, y = expression(theta))
theta.map.1Mb <- theta.map.1Mb + theme(text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

maps.1Mb <- plot_grid(rho.map.1Mb, theta.map.1Mb, nrow = 2, ncol = 1)

# Real data 2L

r2.dm.tab <- as.data.frame(matrix(ncol = 6, nrow = 3))
colnames(r2.dm.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA", "Bin_size(kb)")

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
# filters based on missing data
dm.lands.50kb <- dm.lands.50kb[which(intersect.50kb == F),]

# OLS

# centering
dm.lands.50kb$thetaC <- dm.lands.50kb$theta - mean(dm.lands.50kb$theta)
dm.lands.50kb$tmrcaC <- dm.lands.50kb$tmrca - mean(dm.lands.50kb$tmrca)
dm.lands.50kb$rhoC <- dm.lands.50kb$rho - mean(dm.lands.50kb$rho)

m.diversity <- lm(diversity ~ thetaC + rhoC + tmrcaC + tmrcaC*thetaC, data = dm.lands.50kb)

summary(m.diversity)

plot(resid(m.diversity)~fitted(m.diversity))
dwtest(m.diversity)
hmctest(m.diversity)
hist(resid(m.diversity))

# type 2 ANOVA
anova.diversity <- Anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.dm.tab[1,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                    anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, anova.diversity$VarExp[4] * 100, 50)

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
# filters based on missing data 
dm.lands.200kb <- dm.lands.200kb[which(intersect.200kb == F),]

# OLS

dm.lands.200kb$thetaC <- dm.lands.200kb$theta - mean(dm.lands.200kb$theta)
dm.lands.200kb$tmrcaC <- dm.lands.200kb$tmrca - mean(dm.lands.200kb$tmrca)
dm.lands.200kb$rhoC <- dm.lands.200kb$rho - mean(dm.lands.200kb$rho)

m.diversity <- lm(diversity ~ thetaC + rhoC + tmrcaC + tmrcaC*thetaC, data = dm.lands.200kb)

summary(m.diversity)

# type 2 ANOVA
anova.diversity <- Anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.dm.tab[2,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                   anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, anova.diversity$VarExp[4] * 100, 200)

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
# filters based on missing data
dm.lands.1Mb <- dm.lands.1Mb[which(intersect.1Mb == F),]

# OLS
dm.lands.1Mb$thetaC <- dm.lands.1Mb$theta - mean(dm.lands.1Mb$theta)
dm.lands.1Mb$tmrcaC <- dm.lands.1Mb$tmrca - mean(dm.lands.1Mb$tmrca)
dm.lands.1Mb$rhoC <- dm.lands.1Mb$rho - mean(dm.lands.1Mb$rho)

m.diversity <- lm(diversity ~ thetaC + rhoC + tmrcaC + tmrcaC*thetaC, data = dm.lands.1Mb)


# type 2 ANOVA
anova.diversity <- Anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)


r2.dm.tab[3,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                   anova.diversity$VarExp[1] * 100, 0.0, anova.diversity$VarExp[3] * 100, anova.diversity$VarExp[4] * 100, 1000)

summary(m.diversity)

# R2 Plot

r2.tab.2 <- as.data.frame(cbind(apply(r2.dm.tab, 2, as.numeric)))
names(r2.tab.2)[6] <- "bin.size"

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
r2.plot <- r2.plot + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16))

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

dm.lands.50kb.2L$chr <- "2L"

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

# filters based on missing data 
dm.lands.50kb.2R <- dm.lands.50kb.2R[which(intersect.50kb.2R == F),]

dm.lands.50kb.2R$chr <- "2R"

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

dm.lands.50kb.3L$chr <- "3L"

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

# filters based on missing data
dm.lands.50kb.3R <- dm.lands.50kb.3R[which(intersect.50kb.3R == F),]

dm.lands.50kb.3R$chr <- "3R"

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
dm.lands.50kb$thetaC <- dm.lands.50kb$theta - mean(dm.lands.50kb$theta)
dm.lands.50kb$tmrcaC <- dm.lands.50kb$tmrca - mean(dm.lands.50kb$tmrca)
dm.lands.50kb$rhoC <- dm.lands.50kb$rho - mean(dm.lands.50kb$rho)

dm.lands.50kb$bin <- 1:nrow(dm.lands.50kb)

m.diversity <- lm(diversity ~ (thetaC + rhoC + tmrcaC + thetaC*tmrcaC) * as.factor(chr), data = dm.lands.50kb)
m.diversity2 <- lm(diversity ~ (thetaC + rhoC)* as.factor(chr), data = dm.lands.50kb)
summary(m.diversity2)

hist(resid(m.diversity), nclass = 30)
hmctest(m.diversity, nsim = 3000) 
Box.test(resid(m.diversity)[order(predict(m.diversity))], type = "Ljung-Box") # ***

plot(diversity ~ tmrca, data = dm.lands.50kb)   
summary(m.diversity)
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
anova.diversity <- Anova(m.diversity)
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
d <- dm.lands.50kb
d$y <- (d$diversity^l - 1) / l
g.diversity <- gls(y ~ (theta + rho + tmrca) * as.factor(chr), data = d, corr = corAR1(0, ~bin | chr))

plot(y=resid(g.diversity, type = "pearson"), x=fitted(g.diversity))

summary(g.diversity)
# Parameter estimate(s):
# Phi m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)

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
# Divergence
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
# Evolutionary (Protein) Rates
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

# back to chr idx
dm.lands.evolrate$chr <- as.character(dm.lands.evolrate$chr)
dm.lands.evolrate$chr[which(dm.lands.evolrate$chr == "2L")] <- 2
dm.lands.evolrate$chr[which(dm.lands.evolrate$chr == "2R")] <- 3
dm.lands.evolrate$chr[which(dm.lands.evolrate$chr == "3L")] <- 4
dm.lands.evolrate$chr[which(dm.lands.evolrate$chr == "3R")] <- 5

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
# True (simulated) landscapes
#
#####################

# 50 kb
sim.rho.50kb <- read.table("Misc/dm.sim.rho.50000.txt", header = T)
sim.theta.50kb <- read.table("Misc/dm.sim.theta.50000.txt", header = T)

# rep 1
rep_1.pi.50kb <- read.table("rep_1/rep_1.pi.50kb.bedgraph", header = T)
rep_1.tmrca.50kb <- read.table("rep_1/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_1.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_1.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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

r2.sim.50kb[1, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 1] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 1] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 1] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 1] <- anova.diversity$VarExp[4] * 100

g.div.50kb.1 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                  data = sim.lands.50kb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

g.div.50kb.2 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                    data = sim.lands.50kb, weights = varPower(0, ~theta), cor = corAR1(0, ~bin), method = "ML")

g.div.50kb.3 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                    data = sim.lands.50kb, weights = varPower(0, ~theta), method = "ML")

g.div.50kb.4 <- gls(diversity ~ thetaC + tmrcaC + rhoC + thetaC:tmrcaC,
                    data = sim.lands.50kb, weights = varPower(0, ~tmrcaC), method = "ML")

AIC(g.div.50kb.1, g.div.50kb.2, g.div.50kb.3, g.div.50kb.4)

summary(g.div.50kb.3)
# Value  Std.Error   t-value p-value
# (Intercept)   0.0206901 0.00001905 1086.3241  0.0000
# thetaC        1.3155987 0.00236880  555.3850  0.0000
# tmrcaC        0.0234411 0.00026760   87.5965  0.0000
# rhoC          0.0087173 0.00554727    1.5715  0.1166
# thetaC:tmrcaC 1.4652196 0.03311231   44.2500  0.0000

vif(g.div.50kb.3)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.003575      1.323896      1.004633      1.325415 

# rep_2
rep_2.pi.50kb <- read.table("rep_2/rep_2.pi.50kb.bedgraph", header = T)
rep_2.tmrca.50kb <- read.table("rep_2/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_2.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_2.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    2.057e-02  2.047e-05 1005.301   <2e-16 ***
# thetaC         1.299e+00  2.448e-03  530.815   <2e-16 ***
# tmrcaC         2.385e-02  2.831e-04   84.273   <2e-16 ***
# rhoC          -1.009e-02  6.919e-03   -1.459    0.145    
# thetaC:tmrcaC  1.515e+00  3.222e-02   47.015   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.009894      1.015137      1.005565      1.018886 

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.50kb[1, 2] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 2] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 2] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 2] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 2] <- anova.diversity$VarExp[4] * 100

# rep_3
rep_3.pi.50kb <- read.table("rep_3/rep_3.pi.50kb.bedgraph", header = T)
rep_3.tmrca.50kb <- read.table("rep_3/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_3.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_3.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.0206720  0.0000187 1105.47   <2e-16 ***
# thetaC        1.3103632  0.0022312  587.29   <2e-16 ***
# tmrcaC        0.0236099  0.0002297  102.80   <2e-16 ***
# rhoC          0.0064418  0.0063164    1.02    0.308    
# thetaC:tmrcaC 1.5012872  0.0254605   58.97   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.005320      1.012606      1.004153      1.011333 

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.50kb[1, 3] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 3] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 3] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 3] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 3] <- anova.diversity$VarExp[4] * 100

# rep_4
rep_4.pi.50kb <- read.table("rep_4/rep_4.pi.50kb.bedgraph", header = T)
rep_4.tmrca.50kb <- read.table("rep_4/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_4.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_4.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.0206720  0.0000187 1105.47   <2e-16 ***
# thetaC        1.3103632  0.0022312  587.29   <2e-16 ***
# tmrcaC        0.0236099  0.0002297  102.80   <2e-16 ***
# rhoC          0.0064418  0.0063164    1.02    0.308    
# thetaC:tmrcaC 1.5012872  0.0254605   58.97   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.002369      1.001952      1.002792      1.001253 

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.50kb[1, 4] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 4] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 4] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 4] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 4] <- anova.diversity$VarExp[4] * 100

# rep_5
rep_5.pi.50kb <- read.table("rep_5/rep_5.pi.50kb.bedgraph", header = T)
rep_5.tmrca.50kb <- read.table("rep_5/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_5.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_5.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.062e-02  1.967e-05 1048.31   <2e-16 ***
# thetaC        1.308e+00  2.348e-03  557.06   <2e-16 ***
# tmrcaC        2.374e-02  2.674e-04   88.80   <2e-16 ***
# rhoC          1.528e-03  6.645e-03    0.23    0.818    
# thetaC:tmrcaC 1.471e+00  3.126e-02   47.05   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.005626      1.003205      1.004046      1.005141 

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.50kb[1, 5] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 5] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 5] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 5] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 5] <- anova.diversity$VarExp[4] * 100

# rep_6
rep_6.pi.50kb <- read.table("rep_6/rep_6.pi.50kb.bedgraph", header = T)
rep_6.tmrca.50kb <- read.table("rep_6/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_6.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_6.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.058e-02  2.073e-05 993.114   <2e-16 ***
# thetaC         1.304e+00  2.479e-03 526.030   <2e-16 ***
# tmrcaC         2.363e-02  2.752e-04  85.861   <2e-16 ***
# rhoC          -1.179e-02  7.005e-03  -1.683   0.0929 .  
# thetaC:tmrcaC  1.439e+00  3.174e-02  45.345   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.008799      1.002153      1.004015      1.008954 

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.50kb[1, 6] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 6] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 6] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 6] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 6] <- anova.diversity$VarExp[4] * 100

# rep_7
rep_7.pi.50kb <- read.table("rep_7/rep_7.pi.50kb.bedgraph", header = T)
rep_7.tmrca.50kb <- read.table("rep_7/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_7.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_7.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    2.061e-02  1.945e-05 1059.614   <2e-16 ***
# thetaC         1.313e+00  2.320e-03  565.892   <2e-16 ***
# tmrcaC         2.389e-02  2.547e-04   93.801   <2e-16 ***
# rhoC          -5.753e-03  6.560e-03   -0.877    0.381    
# thetaC:tmrcaC  1.491e+00  3.019e-02   49.398   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC
# 1.005903      1.005336      1.002241      1.003949 

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.50kb[1, 7] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 7] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 7] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 7] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 7] <- anova.diversity$VarExp[4] * 100


# rep_8
rep_8.pi.50kb <- read.table("rep_8/rep_8.pi.50kb.bedgraph", header = T)
rep_8.tmrca.50kb <- read.table("rep_8/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_8.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_8.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.057e-02  2.065e-05 996.284   <2e-16 ***
# thetaC         1.303e+00  2.468e-03 527.803   <2e-16 ***
# tmrcaC         2.380e-02  2.848e-04  83.583   <2e-16 ***
# rhoC          -8.730e-04  6.970e-03  -0.125      0.9    
# thetaC:tmrcaC  1.383e+00  3.162e-02  43.722   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.010112      1.021487      1.003974      1.023819 

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.50kb[1, 8] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 8] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 8] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 8] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 8] <- anova.diversity$VarExp[4] * 100

# rep_9
rep_9.pi.50kb <- read.table("rep_9/rep_9.pi.50kb.bedgraph", header = T)
rep_9.tmrca.50kb <- read.table("rep_9/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_9.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_9.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

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
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)    2.065e-02  1.899e-05 1087.271   <2e-16 ***
# thetaC         1.308e+00  2.265e-03  577.567   <2e-16 ***
# tmrcaC         2.404e-02  2.565e-04   93.731   <2e-16 ***
# rhoC          -7.769e-03  6.412e-03   -1.212    0.226    
# thetaC:tmrcaC  1.508e+00  2.888e-02   52.235   <2e-16 ***

vif(m.div.50kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.003519      1.001713      1.002537      1.001897 

anova.diversity <- Anova(m.div.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.50kb[1, 10] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 10] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 10] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 10] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 10] <- anova.diversity$VarExp[4] * 100



