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

# R_2 table for plotting at the end (real drosophila data)
r2.dm.tab <- as.data.frame(matrix(ncol = 6, nrow = 3))
colnames(r2.dm.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA", "bin.size")

theme.blank <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nreps <- 10

reps <- character(length = nreps)
for(i in 1:nreps) {
  reps[i] <- paste("rep_", i, sep = "")
}

##################c####################
#
# Drosophila-like neutral simulations of 2L (Inferred landscapes)
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


g.rep_1.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep1, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_1.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept) 0.0206321 0.00021240 97.13730  0.0000
# thetaC      1.1816867 0.01703572 69.36522  0.0000
# rhoC        0.0057369 0.08137302  0.07050  0.9438

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

g.rep_2.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep2, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_2.no.tmrca)
# Value  Std.Error   t-value p-value
# (Intercept)  0.0205980 0.00016691 123.40674  0.0000
# thetaC       1.1278261 0.01545460  72.97670  0.0000
# rhoC        -0.0736804 0.08844426  -0.83307  0.4051

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


g.rep_3.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep3, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_3.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0207420 0.00021384 96.99954  0.0000
# thetaC       1.1398975 0.01768070 64.47130  0.0000
# rhoC        -0.1017491 0.09007408 -1.12962  0.2591

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

g.rep_4.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep4, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_4.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0206892 0.00021819 94.82308  0.0000
# thetaC       1.1585830 0.01825400 63.47010  0.0000
# rhoC        -0.1715986 0.09086145 -1.88857  0.0594

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


g.rep_5.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep5, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_5.no.tmrca)
# Value  Std.Error   t-value p-value
# (Intercept)  0.0205709 0.00019283 106.68153  0.0000
# thetaC       1.1335447 0.01727411  65.62102  0.0000
# rhoC        -0.1554325 0.09617437  -1.61615  0.1066

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

g.rep_6.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep6, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_6.no.tmrca)
# Value  Std.Error   t-value p-value
# (Intercept)  0.0205846 0.00018261 112.72154  0.0000
# thetaC       1.1129166 0.01585204  70.20654  0.0000
# rhoC        -0.1062152 0.08342243  -1.27322  0.2034


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

g.rep_7.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep7, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_7.no.tmrca)
# Value  Std.Error   t-value p-value
# (Intercept) 0.0205801 0.00019340 106.41030  0.0000
# thetaC      1.1283864 0.01620181  69.64571  0.0000
# rhoC        0.0166423 0.08056264   0.20658  0.8364

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


g.rep_8.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep8, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_8.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0208143 0.00020979 99.21612  0.0000
# thetaC       1.1333322 0.01746884 64.87738  0.0000
# rhoC        -0.0900602 0.08770082 -1.02690  0.3049

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

g.rep_9.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.50k.rep9, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_9.no.tmrca)
# Value  Std.Error   t-value p-value
# (Intercept)  0.0205299 0.00017580 116.77810  0.0000
# thetaC       1.1341848 0.01570929  72.19836  0.0000
# rhoC        -0.1665026 0.09061323  -1.83751  0.0666

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


g.rep_10.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                         data = inf.lands.50k.rep10, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_10.no.tmrca)
# Value  Std.Error   t-value p-value
# (Intercept)  0.0207082 0.00019363 106.94544  0.0000
# thetaC       1.1606422 0.01631116  71.15634  0.0000
# rhoC        -0.1660896 0.08076784  -2.05638  0.0402

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


g.rep_1.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep1, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_1.no.tmrca)
# Value  Std.Error   t-value p-value
# (Intercept)  0.0206297 0.00020572 100.27847  0.0000
# thetaC       1.1244227 0.02111497  53.25239  0.0000
# rhoC        -0.5711032 0.21962397  -2.60037  0.0103


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

g.rep_2.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep2, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_2.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0206008 0.00021211 97.12366  0.0000
# thetaC       1.1262748 0.01942595 57.97785  0.0000
# rhoC        -0.3984795 0.20018010 -1.99061  0.0484

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


g.rep_3.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep3, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_3.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0207429 0.00024772 83.73477   0.000
# thetaC       1.1273552 0.02313319 48.73324   0.000
# rhoC        -0.5122430 0.25571026 -2.00322   0.047

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

g.rep_4.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep4, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_4.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0206898 0.00024597 84.11461  0.0000
# thetaC       1.1309856 0.02340816 48.31587  0.0000
# rhoC        -0.4082832 0.24994231 -1.63351  0.1045

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

g.rep_5.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep5, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_5.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0205716 0.00023791 86.46709  0.0000
# thetaC       1.1258771 0.02302967 48.88811  0.0000
# rhoC        -0.6236980 0.25913685 -2.40683  0.0173

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

g.rep_6.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep6, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_6.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0205840 0.00022012 93.51378   0e+00
# thetaC       1.1123598 0.02006909 55.42651   0e+00
# rhoC        -0.7401946 0.19645270 -3.76780   2e-04

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


g.rep_7.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep7, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_7.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0205796 0.00021783 94.47441  0.0000
# thetaC       1.1122060 0.02050319 54.24550  0.0000
# rhoC        -0.4602880 0.22027764 -2.08958  0.0384

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


g.rep_8.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep8, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_8.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0208111 0.00023518 88.49129  0.0000
# thetaC       1.1297636 0.02165099 52.18069  0.0000
# rhoC        -0.4075453 0.22773180 -1.78958  0.0756

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


g.rep_9.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.200k.rep9, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_9.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0205379 0.00022295 92.11977  0.0000
# thetaC       1.1066746 0.01936981 57.13401  0.0000
# rhoC        -0.4865608 0.21452090 -2.26813  0.0248

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


g.rep_10.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                         data = inf.lands.200k.rep10, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_10.no.tmrca)
# Value  Std.Error  t-value p-value
# (Intercept)  0.0207131 0.00023459 88.29634  0.0000
# thetaC       1.1352063 0.02089711 54.32360  0.0000
# rhoC        -0.5008572 0.21144100 -2.36878  0.0191

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

g.rep_1.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep1, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_1.no.tmrca)
# Value Std.Error   t-value p-value
# (Intercept)  0.0206090 0.0001242 165.95945  0.0000
# thetaC       1.0935838 0.0263685  41.47310  0.0000
# rhoC        -0.9523384 0.3871228  -2.46004  0.0206

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


g.rep_2.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep2, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_2.no.tmrca)
# Value Std.Error   t-value p-value
# (Intercept)  0.0205937 0.0001941 106.12154  0.0000
# thetaC       1.0716533 0.0343149  31.23000  0.0000
# rhoC        -0.5592219 0.6546032  -0.85429  0.4005

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


g.rep_3.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep3, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_3.no.tmrca)
# Value Std.Error  t-value p-value
# (Intercept)  0.0207330 0.0002082 99.56794  0.0000
# thetaC       1.1012864 0.0377084 29.20532  0.0000
# rhoC        -0.3177321 0.6160273 -0.51578  0.6102

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

g.rep_4.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep4, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_4.no.tmrca)
# Value Std.Error   t-value p-value
# (Intercept)  0.0206798 0.0001348 153.42630  0.0000
# thetaC       1.1175389 0.0287817  38.82807  0.0000
# rhoC        -1.1557791 0.4865494  -2.37546  0.0249

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

g.rep_5.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep5, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_5.no.tmrca)
# Value Std.Error   t-value p-value
# (Intercept)  0.0205697 0.0001476 139.31706  0.0000
# thetaC       1.0982569 0.0293493  37.42021  0.0000
# rhoC        -1.2150813 0.4706686  -2.58161  0.0156

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


g.rep_6.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep6, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_6.no.tmrca)
# Value Std.Error   t-value p-value
# (Intercept)  0.0205747 0.0001889 108.91122  0.0000
# thetaC       1.0881740 0.0343243  31.70273  0.0000
# rhoC        -0.8461274 0.4782870  -1.76908  0.0882

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

g.rep_7.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep7, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_7.no.tmrca)
# Value Std.Error  t-value p-value
# (Intercept)  0.0205818 0.0002149 95.79213  0.0000
# thetaC       1.0897560 0.0336214 32.41258  0.0000
# rhoC        -0.8110199 0.5672255 -1.42980  0.1642

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

g.rep_8.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep8, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_8.no.tmrca)
# Value Std.Error   t-value p-value
# (Intercept)  0.0208035 0.0001699 122.42339  0.0000
# thetaC       1.0823875 0.0328438  32.95562  0.0000
# rhoC        -0.5688595 0.5260339  -1.08141  0.2891

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

g.rep_9.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                        data = inf.lands.1M.rep9, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_9.no.tmrca)
# Value Std.Error   t-value p-value
# (Intercept)  0.0205144 0.0001444 142.03128  0.0000
# thetaC       1.0817758 0.0289649  37.34776  0.0000
# rhoC        -0.9433213 0.4654231  -2.02680  0.0527


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

g.rep_10.no.tmrca <- gls(diversity ~ thetaC + rhoC,
                         data = inf.lands.1M.rep10, cor = corAR1(0, ~bin), method = "ML")

summary(g.rep_10.no.tmrca)
# Value Std.Error   t-value p-value
# (Intercept)  0.020705 0.0001980 104.59296  0.0000
# thetaC       1.062802 0.0337999  31.44392  0.0000
# rhoC        -0.727976 0.5917022  -1.23031  0.2292


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

######################################
#
# Real Drosophila data 50kb windows
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

dm.lands.50kb.2L$thetaC <- dm.lands.50kb.2L$theta - mean(dm.lands.50kb.2L$theta)
dm.lands.50kb.2L$tmrcaC <- dm.lands.50kb.2L$tmrca - mean(dm.lands.50kb.2L$tmrca)
dm.lands.50kb.2L$rhoC <- dm.lands.50kb.2L$rho - mean(dm.lands.50kb.2L$rho)

g.div.dm.50kb.2L <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.50kb.2L, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.50kb.2L)
# Value  Std.Error  t-value p-value
# (Intercept)   0.0097130 0.00001087 893.9585  0.0000
# thetaC        0.9874554 0.00504783 195.6197  0.0000
# rhoC          0.0017748 0.00146095   1.2148  0.2253
# tmrcaC        0.0126819 0.00021188  59.8534  0.0000
# thetaC:tmrcaC 1.2080584 0.04507126  26.8033  0.0000

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

dm.lands.50kb.2R$thetaC <- dm.lands.50kb.2R$theta- mean(dm.lands.50kb.2R$theta)
dm.lands.50kb.2R$tmrcaC <- dm.lands.50kb.2R$tmrca - mean(dm.lands.50kb.2R$tmrca)
dm.lands.50kb.2R$rhoC <- dm.lands.50kb.2R$rho - mean(dm.lands.50kb.2R$rho)

g.div.dm.50kb.2R <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.50kb.2R, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.50kb.2R)
# Value  Std.Error  t-value p-value
# (Intercept)   0.0085652 0.00000938 913.0793  0.0000
# thetaC        0.9709055 0.00361092 268.8801  0.0000
# rhoC          0.0001728 0.00149052   0.1160  0.9078
# tmrcaC        0.0116837 0.00019594  59.6303  0.0000
# thetaC:tmrcaC 1.0680644 0.04933969  21.6472  0.0000

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

dm.lands.50kb.3L$thetaC <- dm.lands.50kb.3L$theta- mean(dm.lands.50kb.3L$theta)
dm.lands.50kb.3L$tmrcaC <- dm.lands.50kb.3L$tmrca - mean(dm.lands.50kb.3L$tmrca)
dm.lands.50kb.3L$rhoC <- dm.lands.50kb.3L$rho - mean(dm.lands.50kb.3L$rho)

g.div.dm.50kb.3L <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.50kb.3L, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.50kb.3L)
# Value   Std.Error  t-value p-value
# (Intercept)   0.0089556 0.000012514 715.6471  0.0000
# thetaC        0.9656448 0.004479703 215.5600  0.0000
# rhoC          0.0028522 0.001435660   1.9867  0.0477
# tmrcaC        0.0118244 0.000153683  76.9400  0.0000
# thetaC:tmrcaC 1.0917168 0.030139577  36.2220  0.0000

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

dm.lands.50kb.3R$thetaC <- dm.lands.50kb.3R$theta- mean(dm.lands.50kb.3R$theta)
dm.lands.50kb.3R$tmrcaC <- dm.lands.50kb.3R$tmrca - mean(dm.lands.50kb.3R$tmrca)
dm.lands.50kb.3R$rhoC <- dm.lands.50kb.3R$rho - mean(dm.lands.50kb.3R$rho)

g.div.dm.50kb.3R <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.50kb.3R, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.50kb.3R)
# Value  Std.Error  t-value p-value
# (Intercept)    0.0068266 0.00001072 636.5458  0.0000
# thetaC         0.9529991 0.00406945 234.1840  0.0000
# rhoC          -0.0020686 0.00150933  -1.3705  0.1715
# tmrcaC         0.0102749 0.00017107  60.0611  0.0000
# thetaC:tmrcaC  1.0361235 0.05282692  19.6136  0.0000

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
# centering
dm.lands.50kb$thetaC <- dm.lands.50kb$theta - mean(dm.lands.50kb$theta)
dm.lands.50kb$tmrcaC <- dm.lands.50kb$tmrca - mean(dm.lands.50kb$tmrca)
dm.lands.50kb$rhoC <- dm.lands.50kb$rho - mean(dm.lands.50kb$rho)

dm.lands.50kb$bin <- 1:nrow(dm.lands.50kb)

m.div.dm.50kb <- lm(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC), data = dm.lands.50kb)

plot(resid(m.div.dm.50kb)~fitted(m.div.dm.50kb))
hist(resid(m.div.dm.50kb), nclass = 30)
dwtest(m.div.dm.50kb)
hmctest(m.div.dm.50kb, nsim = 10000) 

summary(m.div.dm.50kb)
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)   8.521e-03  5.398e-06 1578.575   <2e-16 ***
# thetaC        9.746e-01  1.983e-03  491.576   <2e-16 ***
# rhoC          1.752e-03  7.737e-04    2.264   0.0237 *  
# tmrcaC        1.139e-02  8.435e-05  135.042   <2e-16 ***
# thetaC:tmrcaC 1.051e+00  1.799e-02   58.418   <2e-16 ***

# type 2 ANOVA
anova.diversity <- Anova(m.div.dm.50kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.dm.tab[1,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                    anova.diversity$VarExp[1] * 100,
                    anova.diversity$VarExp[2] * 100,
                    anova.diversity$VarExp[3] * 100,
                    anova.diversity$VarExp[4] * 100,
                    50)

# GLS
g.div.dm.50kb.1 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.50kb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

g.div.dm.50kb.2 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.50kb, weights = varPower(0, ~theta), cor = corAR1(0, ~bin), method = "ML")

g.div.dm.50kb.3 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.50kb, weights = varPower(0, ~theta), method = "ML")

g.div.dm.50kb.4 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.50kb, weights = varPower(0, ~tmrcaC), method = "ML")

AIC(g.div.dm.50kb.1, g.div.dm.50kb.2, g.div.dm.50kb.3, g.div.dm.50kb.4)

summary(g.div.dm.50kb.1)
# Value   Std.Error   t-value p-value
# (Intercept)   0.0085133 0.000006325 1345.9050  0.0000
# thetaC        0.9730593 0.002301256  422.8384  0.0000
# rhoC          0.0014472 0.000798961    1.8114  0.0703
# tmrcaC        0.0114980 0.000091587  125.5412  0.0000
# thetaC:tmrcaC 1.0777123 0.019925270   54.0877  0.0000


# Linear model without TMRCA --> rho becomes significant
g.div.dm.50kb.5 <- gls(diversity ~ (thetaC + rhoC),
                       data = dm.lands.50kb, weights = varPower(0, ~thetaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.50kb.5)
# Value   Std.Error  t-value p-value
# (Intercept) 0.0086582 0.000023593 366.9880       0
# thetaC      1.0990974 0.008957997 122.6945       0
# rhoC        0.0528973 0.002560827  20.6563       0

######################################
#
# Real Drosophila data 200kb windows
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

dm.lands.200kb.2L$chr <- "2L"

dm.lands.200kb.2L$thetaC <- dm.lands.200kb.2L$theta - mean(dm.lands.200kb.2L$theta)
dm.lands.200kb.2L$tmrcaC <- dm.lands.200kb.2L$tmrca - mean(dm.lands.200kb.2L$tmrca)
dm.lands.200kb.2L$rhoC <- dm.lands.200kb.2L$rho - mean(dm.lands.200kb.2L$rho)

g.div.dm.200kb.2L <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.200kb.2L, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.200kb.2L)
# Value  Std.Error  t-value p-value
# (Intercept)   0.0096851 0.00001590 609.2540  0.0000
# thetaC        0.9912075 0.00852140 116.3197  0.0000
# rhoC          0.0035710 0.00328708   1.0864  0.2807
# tmrcaC        0.0122532 0.00061930  19.7855  0.0000
# thetaC:tmrcaC 1.0357572 0.10202136  10.1524  0.0000


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

# filters based on missing data 
dm.lands.200kb.2R <- dm.lands.200kb.2R[which(intersect.200kb.2R == F),]

dm.lands.200kb.2R$chr <- "2R"

dm.lands.200kb.2R$thetaC <- dm.lands.200kb.2R$theta- mean(dm.lands.200kb.2R$theta)
dm.lands.200kb.2R$tmrcaC <- dm.lands.200kb.2R$tmrca - mean(dm.lands.200kb.2R$tmrca)
dm.lands.200kb.2R$rhoC <- dm.lands.200kb.2R$rho - mean(dm.lands.200kb.2R$rho)

g.div.dm.200kb.2R <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.200kb.2R, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.200kb.2R)
# Value  Std.Error  t-value p-value
# (Intercept)   0.0085392 0.00001264 675.8307  0.0000
# thetaC        0.9699152 0.00505288 191.9530  0.0000
# rhoC          0.0027465 0.00292834   0.9379  0.3516
# tmrcaC        0.0110299 0.00043604  25.2959  0.0000
# thetaC:tmrcaC 0.9373743 0.11041693   8.4894  0.0000

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

dm.lands.200kb.3L$chr <- "3L"

dm.lands.200kb.3L$thetaC <- dm.lands.200kb.3L$theta- mean(dm.lands.200kb.3L$theta)
dm.lands.200kb.3L$tmrcaC <- dm.lands.200kb.3L$tmrca - mean(dm.lands.200kb.3L$tmrca)
dm.lands.200kb.3L$rhoC <- dm.lands.200kb.3L$rho - mean(dm.lands.200kb.3L$rho)

g.div.dm.200kb.3L <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.200kb.3L, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.200kb.3L)
# Value  Std.Error  t-value p-value
# (Intercept)   0.0089474 0.00001816 492.7438  0.0000
# thetaC        0.9811147 0.00736731 133.1713  0.0000
# rhoC          0.0078270 0.00261567   2.9923  0.0037
# tmrcaC        0.0106115 0.00042145  25.1784  0.0000
# thetaC:tmrcaC 0.7871946 0.07703046  10.2193  0.0000

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

# filters based on missing data
dm.lands.200kb.3R <- dm.lands.200kb.3R[which(intersect.200kb.3R == F),]

dm.lands.200kb.3R$chr <- "3R"

dm.lands.200kb.3R$thetaC <- dm.lands.200kb.3R$theta- mean(dm.lands.200kb.3R$theta)
dm.lands.200kb.3R$tmrcaC <- dm.lands.200kb.3R$tmrca - mean(dm.lands.200kb.3R$tmrca)
dm.lands.200kb.3R$rhoC <- dm.lands.200kb.3R$rho - mean(dm.lands.200kb.3R$rho)

g.div.dm.200kb.3R <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.200kb.3R, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.200kb.3R)
# (Intercept)   0.0069082 0.00001124 614.5051  0.0000
# thetaC        0.9656153 0.00430329 224.3901  0.0000
# rhoC          0.0012400 0.00220961   0.5612  0.5763
# tmrcaC        0.0091446 0.00030990  29.5085  0.0000
# thetaC:tmrcaC 0.5211798 0.09505901   5.4827  0.0000

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
# 0.20 p-value = 2e-13
cor.test(~rho + tmrca, data = dm.lands.200kb, method = "spearman")
# 0.48 p-value < 2.2e-16
cor.test(~theta + tmrca, data = dm.lands.200kb, method = "spearman")
# 0.46 p-value < 2.2e-16

# Linear models
# centering
dm.lands.200kb$thetaC <- dm.lands.200kb$theta - mean(dm.lands.200kb$theta)
dm.lands.200kb$tmrcaC <- dm.lands.200kb$tmrca - mean(dm.lands.200kb$tmrca)
dm.lands.200kb$rhoC <- dm.lands.200kb$rho - mean(dm.lands.200kb$rho)

dm.lands.200kb$bin <- 1:nrow(dm.lands.200kb)

m.div.dm.200kb <- lm(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC), data = dm.lands.200kb)

plot(resid(m.div.dm.200kb)~fitted(m.div.dm.200kb))
hist(resid(m.div.dm.200kb), nclass = 30)
dwtest(m.div.dm.200kb)
hmctest(m.div.dm.200kb, nsim = 10000) 

summary(m.div.dm.200kb)
# (Intercept)   8.501e-03  8.291e-06 1025.258   <2e-16 ***
# thetaC        9.767e-01  3.253e-03  300.252   <2e-16 ***
# rhoC          3.605e-03  1.546e-03    2.331   0.0204 *  
# tmrcaC        1.113e-02  2.219e-04   50.141   <2e-16 ***
# thetaC:tmrcaC 8.880e-01  4.677e-02   18.987   <2e-16 ***
  
# type 2 ANOVA
anova.diversity <- Anova(m.div.dm.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.dm.tab[2,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                   anova.diversity$VarExp[1] * 100,
                   anova.diversity$VarExp[2] * 100,
                   anova.diversity$VarExp[3] * 100,
                   anova.diversity$VarExp[4] * 100,
                   200)

# GLS
g.div.dm.200kb.1 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.200kb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

g.div.dm.200kb.2 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.200kb, weights = varPower(0, ~thetaC), cor = corAR1(0, ~bin), method = "ML")

g.div.dm.200kb.3 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.200kb, weights = varPower(0, ~thetaC), method = "ML")

g.div.dm.200kb.4 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.200kb, weights = varPower(0, ~tmrcaC), method = "ML")

AIC(g.div.dm.200kb.1, g.div.dm.200kb.2, g.div.dm.200kb.3, g.div.dm.200kb.4)

summary(g.div.dm.200kb.1)
# Value  Std.Error  t-value p-value
# Intercept)   0.0084988 0.00000935 908.9635  0.0000
# thetaC        0.9728966 0.00365565 266.1347  0.0000
# rhoC          0.0036066 0.00165622   2.1776  0.0302
# tmrcaC        0.0110790 0.00023178  47.8001  0.0000
# thetaC:tmrcaC 0.8919768 0.05097523  17.4982  0.0000


# Linear model without TMRCA --> rho becomes significant
g.div.dm.200kb.5 <- gls(diversity ~ (thetaC + rhoC),
                       data = dm.lands.200kb, weights = varPower(0, ~thetaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.200kb.5)
# Value   Std.Error   t-value p-value
# (Intercept) 0.0086072 0.000036015 238.98914       0
# thetaC      1.0796433 0.012533437  86.14104       0
# rhoC        0.0600940 0.004429202  13.56767       0

######################################
#
# Real Drosophila data 1Mb windows
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

# filters based on missing data 
dm.lands.1Mb.2L <- dm.lands.1Mb.2L[which(intersect.1Mb.2L == F),]

dm.lands.1Mb.2L$chr <- "2L"

dm.lands.1Mb.2L$thetaC <- dm.lands.1Mb.2L$theta - mean(dm.lands.1Mb.2L$theta)
dm.lands.1Mb.2L$tmrcaC <- dm.lands.1Mb.2L$tmrca - mean(dm.lands.1Mb.2L$tmrca)
dm.lands.1Mb.2L$rhoC <- dm.lands.1Mb.2L$rho - mean(dm.lands.1Mb.2L$rho)

g.div.dm.1Mb.2L <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.1Mb.2L, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.1Mb.2L)
# Value  Std.Error  t-value p-value
# (Intercept)   0.0101124 0.00002513 402.3734  0.0000
# thetaC        0.9865547 0.02033816  48.5076  0.0000
# rhoC          0.0010558 0.00560117   0.1885  0.8532
# tmrcaC        0.0129654 0.00122549  10.5798  0.0000
# thetaC:tmrcaC 1.8344985 0.31522351   5.8197  0.0000

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

# filters based on missing data 
dm.lands.1Mb.2R <- dm.lands.1Mb.2R[which(intersect.1Mb.2R == F),]

dm.lands.1Mb.2R$chr <- "2R"

dm.lands.1Mb.2R$thetaC <- dm.lands.1Mb.2R$theta- mean(dm.lands.1Mb.2R$theta)
dm.lands.1Mb.2R$tmrcaC <- dm.lands.1Mb.2R$tmrca - mean(dm.lands.1Mb.2R$tmrca)
dm.lands.1Mb.2R$rhoC <- dm.lands.1Mb.2R$rho - mean(dm.lands.1Mb.2R$rho)

g.div.dm.1Mb.2R <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.1Mb.2R, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.1Mb.2R)
# Value Std.Error  t-value p-value
# (Intercept)   0.0091552 0.0000283 323.2402  0.0000
# thetaC        0.9718147 0.0212995  45.6262  0.0000
# rhoC          0.0083025 0.0084990   0.9769  0.3479
# tmrcaC        0.0103641 0.0019099   5.4266  0.0002
# thetaC:tmrcaC 0.5202501 0.5098563   1.0204  0.3277

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

# filters based on missing data ( > 25% per window)
dm.lands.1Mb.3L <- dm.lands.1Mb.3L[which(intersect.1Mb.3L == F),]

dm.lands.1Mb.3L$chr <- "3L"

dm.lands.1Mb.3L$thetaC <- dm.lands.1Mb.3L$theta- mean(dm.lands.1Mb.3L$theta)
dm.lands.1Mb.3L$tmrcaC <- dm.lands.1Mb.3L$tmrca - mean(dm.lands.1Mb.3L$tmrca)
dm.lands.1Mb.3L$rhoC <- dm.lands.1Mb.3L$rho - mean(dm.lands.1Mb.3L$rho)

g.div.dm.1Mb.3L <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.1Mb.3L, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.1Mb.3L)
# Value  Std.Error  t-value p-value
# (Intercept)   0.0091813 0.00002405 381.7940  0.0000
# thetaC        1.0032652 0.01244385  80.6233  0.0000
# rhoC          0.0070903 0.00396199   1.7896  0.0913
# tmrcaC        0.0084717 0.00094627   8.9528  0.0000
# thetaC:tmrcaC 0.3950013 0.16662491   2.3706  0.0298

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

# filters based on missing data
dm.lands.1Mb.3R <- dm.lands.1Mb.3R[which(intersect.1Mb.3R == F),]

dm.lands.1Mb.3R$chr <- "3R"

dm.lands.1Mb.3R$thetaC <- dm.lands.1Mb.3R$theta- mean(dm.lands.1Mb.3R$theta)
dm.lands.1Mb.3R$tmrcaC <- dm.lands.1Mb.3R$tmrca - mean(dm.lands.1Mb.3R$tmrca)
dm.lands.1Mb.3R$rhoC <- dm.lands.1Mb.3R$rho - mean(dm.lands.1Mb.3R$rho)

g.div.dm.1Mb.3R <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                        data = dm.lands.1Mb.3R, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.1Mb.3R)
# Value  Std.Error   t-value p-value
# (Intercept)   0.0067246 0.00002252 298.66683  0.0000
# thetaC        0.9664707 0.00658600 146.74629  0.0000
# rhoC          0.0031144 0.00241054   1.29200  0.2173
# tmrcaC        0.0094832 0.00049711  19.07661  0.0000
# thetaC:tmrcaC 0.7457601 0.17910164   4.16389  0.0010

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
# 0.20 p-value = 2e-13
cor.test(~rho + tmrca, data = dm.lands.1Mb, method = "spearman")
# 0.48 p-value < 2.2e-16
cor.test(~theta + tmrca, data = dm.lands.1Mb, method = "spearman")
# 0.46 p-value < 2.2e-16

# Linear models
# centering
dm.lands.1Mb$thetaC <- dm.lands.1Mb$theta - mean(dm.lands.1Mb$theta)
dm.lands.1Mb$tmrcaC <- dm.lands.1Mb$tmrca - mean(dm.lands.1Mb$tmrca)
dm.lands.1Mb$rhoC <- dm.lands.1Mb$rho - mean(dm.lands.1Mb$rho)

dm.lands.1Mb$bin <- 1:nrow(dm.lands.1Mb)

m.div.dm.1Mb <- lm(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC), data = dm.lands.1Mb)

plot(resid(m.div.dm.1Mb)~fitted(m.div.dm.1Mb))
hist(resid(m.div.dm.1Mb), nclass = 30)
dwtest(m.div.dm.1Mb)
hmctest(m.div.dm.1Mb, nsim = 10000) 

summary(m.div.dm.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   8.808e-03  1.551e-05 568.036  < 2e-16 ***
# thetaC        9.919e-01  6.725e-03 147.497  < 2e-16 ***
# rhoC          6.510e-03  3.034e-03   2.145   0.0353 *  
# tmrcaC        9.525e-03  6.035e-04  15.784  < 2e-16 ***
# thetaC:tmrcaC 4.881e-01  1.131e-01   4.315 4.99e-05 ***

# type 2 ANOVA
anova.diversity <- Anova(m.div.dm.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.dm.tab[3,] <- c((anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100,
                   anova.diversity$VarExp[1] * 100,
                   anova.diversity$VarExp[2] * 100,
                   anova.diversity$VarExp[3] * 100,
                   anova.diversity$VarExp[4] * 100,
                   1000)

# GLS
g.div.dm.1Mb.1 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.1Mb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

g.div.dm.1Mb.2 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.1Mb, weights = varPower(0, ~theta), cor = corAR1(0, ~bin), method = "ML")

g.div.dm.1Mb.3 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.1Mb, weights = varPower(0, ~theta), method = "ML")

g.div.dm.1Mb.4 <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC),
                       data = dm.lands.1Mb, weights = varPower(0, ~tmrcaC), method = "ML")

AIC(g.div.dm.1Mb.1, g.div.dm.1Mb.2, g.div.dm.1Mb.3, g.div.dm.1Mb.4)

summary(g.div.dm.1Mb.1)
# Value   Std.Error   t-value p-value
# (Intercept)   0.0087767 0.00001307 671.7413  0.0000
# thetaC        0.9802560 0.00621327 157.7682  0.0000
# rhoC          0.0011873 0.00289132   0.4106  0.6826
# tmrcaC        0.0104410 0.00061019  17.1110  0.0000
# thetaC:tmrcaC 0.6526937 0.13373906   4.8804  0.0000


# Linear model without TMRCA --> rho becomes significant
g.div.dm.1Mb.5 <- gls(diversity ~ (thetaC + rhoC),
                       data = dm.lands.1Mb, weights = varPower(0, ~thetaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.div.dm.1Mb.5)
# Value   Std.Error  t-value p-value
# (Intercept) 0.0086582 0.000023593 366.9880       0
# thetaC      1.0990974 0.008957997 122.6945       0
# rhoC        0.0528973 0.002560827  20.6563       0

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
          z = lands.divergence.dm$theta, method = "spearman")
# 0.05 p-value = 0.17


########################################
#
# Evolutionary (Protein) Rates
#
########################################

# loads 
dm.raw <- read.table("~/Data/iSMC/theta_paper/real_data/dm_misc/dpgp3_Dyak_bpp.all.csv", header = T, fill = T, stringsAsFactors = T)
dm.tbl <- na.omit(dm.raw)

# gets ratios 
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

write.table(dm.lands.evolrate, "dm_chr_maps/dm_maps_50kb_protein_rates.tsv",
            sep = "\t", quote = F, col.names = T, row.names = F)

# linear model in coding regions
# centering
dm.lands.evolrate$thetaC <- dm.lands.evolrate$theta - mean(dm.lands.evolrate$theta)
dm.lands.evolrate$tmrcaC <- dm.lands.evolrate$tmrca - mean(dm.lands.evolrate$tmrca)
dm.lands.evolrate$rhoC <- dm.lands.evolrate$rho - mean(dm.lands.evolrate$rho)

dm.lands.evolrate$bin <- 1:nrow(dm.lands.evolrate)

m.dm.cds <- lm(diversity ~ (thetaC + rhoC + tmrcaC + thetaC*tmrcaC) * chr, data = dm.lands.evolrate)

plot(resid(m.dm.cds)~fitted(m.dm.cds))
hist(resid(m.dm.cds), nclass = 30)
dwtest(m.dm.cds)
hmctest(m.dm.cds, nsim = 10000) 

summary(m.dm.cds)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          8.597e-03  1.623e-05 529.546  < 2e-16 ***
# thetaC               9.598e-01  6.449e-03 148.817  < 2e-16 ***
# rhoC                 2.041e-03  1.694e-03   1.205 0.228460    
# tmrcaC               1.148e-02  2.104e-04  54.544  < 2e-16 ***
# chr2R               -1.034e-04  1.987e-05  -5.205 2.23e-07 ***
# chr3L               -1.199e-04  1.910e-05  -6.276 4.62e-10 ***
# chr3R               -8.636e-05  2.136e-05  -4.043 5.57e-05 ***
# thetaC:tmrcaC        1.157e+00  4.321e-02  26.772  < 2e-16 ***
# thetaC:chr2R        -9.775e-03  7.982e-03  -1.225 0.220941    
# thetaC:chr3L         1.911e-03  7.237e-03   0.264 0.791752    
# thetaC:chr3R         2.881e-03  7.620e-03   0.378 0.705367    
# rhoC:chr2R          -1.540e-03  2.375e-03  -0.649 0.516759    
# rhoC:chr3L           1.806e-03  2.156e-03   0.838 0.402292    
# rhoC:chr3R          -4.121e-03  2.542e-03  -1.621 0.105232    
# tmrcaC:chr2R         5.242e-04  2.883e-04   1.818 0.069253 .  
# tmrcaC:chr3L        -1.434e-04  2.623e-04  -0.547 0.584648    
# tmrcaC:chr3R         2.060e-04  2.981e-04   0.691 0.489652    
# thetaC:tmrcaC:chr2R -1.827e-01  6.123e-02  -2.985 0.002889 ** 
# thetaC:tmrcaC:chr3L -4.053e-03  5.188e-02  -0.078 0.937754    
# thetaC:tmrcaC:chr3R -2.334e-01  7.070e-02  -3.301 0.000989 ***

# type 2 anova
anova.diversity.cds <- Anova(m.dm.cds)
apiss <- anova.diversity.cds$"Sum Sq"
anova.diversity.cds$VarExp <- apiss / sum(apiss)

anova.diversity.cds
# Sum Sq   Df    F value   Pr(>F)  VarExp
# thetaC            0.0068959    1 2.1893e+05 0.000000 0.92243
# rhoC              0.0000001    1 3.9590e+00 0.046814 0.00002
# tmrcaC            0.0004327    1 1.3739e+04 0.000000 0.05789
# chr               0.0000030    3 3.1741e+01 0.000000 0.00040
# thetaC:tmrcaC     0.0000970    1 3.0806e+03 0.000000 0.01298
# thetaC:chr        0.0000002    3 2.3316e+00 0.072540 0.00003
# rhoC:chr          0.0000002    3 2.3713e+00 0.068825 0.00003
# tmrcaC:chr        0.0000015    3 1.5454e+01 0.000000 0.00020
# thetaC:tmrcaC:chr 0.0000007    3 7.7198e+00 0.000041 0.00010
# Residuals         0.0000444 1409                     0.00594

# GLS
dm.lands.evolrate$bin <- 1:nrow(dm.lands.evolrate)

g.dm.cds <- gls(diversity ~ (thetaC + rhoC + tmrcaC + thetaC:tmrcaC) * chr,
                 data = dm.lands.evolrate, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

summary(g.dm.cds)
# Value  Std.Error  t-value p-value
# (Intercept)          0.0085814 0.00001683 509.9168  0.0000
# thetaC               0.9604655 0.00602716 159.3561  0.0000
# rhoC                 0.0017647 0.00164068   1.0756  0.2823
# tmrcaC               0.0115725 0.00021044  54.9909  0.0000
# chr2R               -0.0000826 0.00002063  -4.0039  0.0001
# chr3L               -0.0001170 0.00001990  -5.8790  0.0000
# chr3R               -0.0000684 0.00002231  -3.0663  0.0022
# thetaC:tmrcaC        1.2054766 0.04727993  25.4966  0.0000
# thetaC:chr2R        -0.0086915 0.00759471  -1.1444  0.2526
# thetaC:chr3L        -0.0081094 0.00679677  -1.1931  0.2330
# thetaC:chr3R         0.0040474 0.00724860   0.5584  0.5767
# rhoC:chr2R          -0.0013860 0.00232361  -0.5965  0.5510
# rhoC:chr3L           0.0006780 0.00208428   0.3253  0.7450
# rhoC:chr3R          -0.0048064 0.00243725  -1.9721  0.0488
# tmrcaC:chr2R         0.0002858 0.00029480   0.9695  0.3325
# tmrcaC:chr3L         0.0000623 0.00025952   0.2400  0.8103
# tmrcaC:chr3R        -0.0000453 0.00031256  -0.1450  0.8847
# thetaC:tmrcaC:chr2R -0.2167039 0.06760625  -3.2054  0.0014
# thetaC:tmrcaC:chr3L -0.0165022 0.05588164  -0.2953  0.7678
# thetaC:tmrcaC:chr3R -0.3127489 0.07676455  -4.0741  0.0000


summary(g.dm.cds.pis)
# correlations
cor.test(dm.lands.evolrate$PiN, dm.lands.evolrate$theta, method = "spearman") 
# 0.08, p-value 0.06
cor.test(dm.lands.evolrate$PiS, dm.lands.evolrate$theta, method = "spearman") 
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

########################################
#
# Drosophila-like neutral simulations of 2L (True landscapes)
#
#######################################

setwd("~")
setwd("Data/iSMC/theta_paper/sim_data/rs_drosophila_2/")

# 50 kb
r2.sim.50kb <- as.data.frame(matrix(ncol = nreps, nrow = 5))
row.names(r2.sim.50kb) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA")
colnames(r2.sim.50kb) <- reps

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

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
m.div.50kb.2 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC + rhoC:tmrcaC, data = sim.lands.50kb)
m.div.50kb.3 <- lm(diversity ~ (thetaC + rhoC + tmrcaC) ^ 2, data = sim.lands.50kb)

AIC(m.div.50kb, m.div.50kb.2, m.div.50kb.3)

plot(resid(m.div.50kb)~fitted(m.div.50kb))
dwtest(m.div.50kb)
hmctest(m.div.50kb)
hist(resid(m.div.50kb))

summary(m.div.50kb)
# Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)   2.069e-02  1.908e-05 1084.381   <2e-16 ***
# thetaC        1.309e+00  2.283e-03  573.279   <2e-16 ***
# tmrcaC        2.342e-02  2.686e-04   87.195   <2e-16 ***
# rhoC          1.435e-02  6.455e-03    2.223   0.0266 *  
# thetaC:tmrcaC 1.483e+00  3.029e-02   48.983   <2e-16 ***

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

g.div.50kb.1 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                  data = sim.lands.50kb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

g.div.50kb.2 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = sim.lands.50kb, weights = varPower(0, ~theta), cor = corAR1(0, ~bin), method = "ML")

g.div.50kb.3 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = sim.lands.50kb, weights = varPower(0, ~theta), method = "ML")

g.div.50kb.4 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
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


cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
m.div.50kb.2 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC + rhoC:tmrcaC, data = sim.lands.50kb)
m.div.50kb.3 <- lm(diversity ~ (thetaC + rhoC + tmrcaC) ^ 2, data = sim.lands.50kb)

AIC(m.div.50kb, m.div.50kb.2, m.div.50kb.3)

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

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
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

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
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

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
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

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
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

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
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

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
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

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
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

r2.sim.50kb[1, 9] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.50kb[2, 9] <- anova.diversity$VarExp[1] * 100
r2.sim.50kb[3, 9] <- anova.diversity$VarExp[2] * 100
r2.sim.50kb[4, 9] <- anova.diversity$VarExp[3] * 100
r2.sim.50kb[5, 9] <- anova.diversity$VarExp[4] * 100


# rep_10
rep_10.pi.50kb <- read.table("rep_10/rep_10.pi.50kb.bedgraph", header = T)
rep_10.tmrca.50kb <- read.table("rep_10/sim.tmrca.50k.map", header = T)

sim.lands.50kb <- as.data.frame(cbind(rep_10.pi.50kb$pi, sim.theta.50kb$sim, sim.rho.50kb$sim, rep_10.tmrca.50kb$tmrca))
names(sim.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.50kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.50kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.50kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.50kb, method = "spearman") 

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

m.div.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.50kb)
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



r2.sim.50kb$average <- rowMeans(r2.sim.50kb)
r2.sim.50kb <- transform(r2.sim.50kb, sd=apply(r2.sim.50kb, 1, sd, na.rm = TRUE))

# 200 kb
r2.sim.200kb <- as.data.frame(matrix(ncol = nreps, nrow = 5))
row.names(r2.sim.200kb) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA")
colnames(r2.sim.200kb) <- reps

sim.rho.200kb <- read.table("Misc/dm.sim.rho.2e+05.txt", header = T)
sim.theta.200kb <- read.table("Misc/dm.sim.theta.2e+05.txt", header = T)

# rep 1
rep_1.pi.200kb <- read.table("rep_1/rep_1.pi.200kb.bedgraph", header = T)
rep_1.tmrca.200kb <- read.table("rep_1/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_1.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_1.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.067e-02  2.788e-05 741.370   <2e-16 ***
# thetaC        1.306e+00  3.769e-03 346.532   <2e-16 ***
# tmrcaC        2.526e-02  7.654e-04  32.997   <2e-16 ***
# rhoC          3.201e-02  1.870e-02   1.712    0.089 .  
# thetaC:tmrcaC 1.536e+00  1.072e-01  14.325   <2e-16 ***

vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.006587      1.000386      1.015701      1.009284 

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 1] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 1] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 1] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 1] <- anova.diversity$VarExp[4] * 100


# rep_2
rep_2.pi.200kb <- read.table("rep_2/rep_2.pi.200kb.bedgraph", header = T)
rep_2.tmrca.200kb <- read.table("rep_2/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_2.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_2.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.056e-02  2.494e-05  824.14   <2e-16 ***
# thetaC         1.297e+00  3.394e-03  382.15   <2e-16 ***
# tmrcaC         2.455e-02  6.157e-04   39.88   <2e-16 ***
# rhoC          -3.118e-02  1.667e-02   -1.87   0.0635 .  
# thetaC:tmrcaC  1.643e+00  8.199e-02   20.04   <2e-16 ***

vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.020463      1.020098      1.009539      1.034557 

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 2] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 2] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 2] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 2] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 2] <- anova.diversity$VarExp[4] * 100

# rep_3
rep_3.pi.200kb <- read.table("rep_3/rep_3.pi.200kb.bedgraph", header = T)
rep_3.tmrca.200kb <- read.table("rep_3/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_3.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_3.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.064e-02  3.245e-05 636.263   <2e-16 ***
# thetaC        1.309e+00  4.394e-03 297.825   <2e-16 ***
# tmrcaC        2.488e-02  7.859e-04  31.662   <2e-16 ***
# rhoC          3.167e-02  2.173e-02   1.458    0.147    
# thetaC:tmrcaC 1.472e+00  1.097e-01  13.422   <2e-16 ***
  
vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.011267      1.007777      1.013794      1.007048

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 3] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 3] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 3] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 3] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 3] <- anova.diversity$VarExp[4] * 100

# rep_4
rep_4.pi.200kb <- read.table("rep_4/rep_4.pi.200kb.bedgraph", header = T)
rep_4.tmrca.200kb <- read.table("rep_4/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_4.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_4.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.0205018  0.0000271 756.559   <2e-16 ***
# thetaC        1.3010756  0.0036626 355.236   <2e-16 ***
# tmrcaC        0.0244436  0.0006814  35.873   <2e-16 ***
# rhoC          0.0087841  0.0181026   0.485    0.628    
# thetaC:tmrcaC 1.4688735  0.0873643  16.813   <2e-16 ***

vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.007460      1.019677      1.009198      1.019709 

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 4] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 4] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 4] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 4] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 4] <- anova.diversity$VarExp[4] * 100

# rep_5
rep_5.pi.200kb <- read.table("rep_5/rep_5.pi.200kb.bedgraph", header = T)
rep_5.tmrca.200kb <- read.table("rep_5/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_5.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_5.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.061e-02  2.367e-05 870.811   <2e-16 ***
# thetaC        1.306e+00  3.205e-03 407.532   <2e-16 ***
# tmrcaC        2.402e-02  6.327e-04  37.957   <2e-16 ***
# rhoC          8.924e-03  1.586e-02   0.563    0.575    
# thetaC:tmrcaC 1.380e+00  7.367e-02  18.728   <2e-16 ***

vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.011705      1.058659      1.015517      1.054433 

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 5] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 5] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 5] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 5] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 5] <- anova.diversity$VarExp[4] * 100

# rep_6
rep_6.pi.200kb <- read.table("rep_6/rep_6.pi.200kb.bedgraph", header = T)
rep_6.tmrca.200kb <- read.table("rep_6/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_6.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_6.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.055e-02  2.337e-05 879.460   <2e-16 ***
# thetaC        1.299e+00  3.288e-03 395.145   <2e-16 ***
# tmrcaC        2.392e-02  6.333e-04  37.771   <2e-16 ***
# rhoC          3.569e-03  1.559e-02   0.229    0.819    
# thetaC:tmrcaC 1.533e+00  6.794e-02  22.563   <2e-16 ***

vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.092293      1.019078      1.006542      1.103193 

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 6] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 6] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 6] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 6] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 6] <- anova.diversity$VarExp[4] * 100

# rep_7
rep_7.pi.200kb <- read.table("rep_7/rep_7.pi.200kb.bedgraph", header = T)
rep_7.tmrca.200kb <- read.table("rep_7/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_7.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_7.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.062e-02  3.141e-05 656.642   <2e-16 ***
# thetaC         1.315e+00  4.235e-03 310.525   <2e-16 ***
# tmrcaC         2.369e-02  8.316e-04  28.486   <2e-16 ***
# rhoC          -1.620e-02  2.076e-02  -0.781    0.436    
# thetaC:tmrcaC  1.422e+00  1.203e-01  11.823   <2e-16 ***

vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.028877      1.029367      1.013152      1.003227 

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 7] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 7] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 7] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 7] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 7] <- anova.diversity$VarExp[4] * 100


# rep_8
rep_8.pi.200kb <- read.table("rep_8/rep_8.pi.200kb.bedgraph", header = T)
rep_8.tmrca.200kb <- read.table("rep_8/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_8.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_8.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.057e-02  2.534e-05 811.819   <2e-16 ***
# thetaC         1.309e+00  3.478e-03 376.496   <2e-16 ***
# tmrcaC         2.531e-02  6.716e-04  37.681   <2e-16 ***
# rhoC          -9.189e-03  1.709e-02  -0.538    0.592    
# thetaC:tmrcaC  1.518e+00  7.926e-02  19.154   <2e-16 ***

vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.046676      1.072450      1.036894      1.089271 

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 8] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 8] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 8] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 8] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 8] <- anova.diversity$VarExp[4] * 100

# rep_9
rep_9.pi.200kb <- read.table("rep_9/rep_9.pi.200kb.bedgraph", header = T)
rep_9.tmrca.200kb <- read.table("rep_9/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_9.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_9.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.061e-02  2.834e-05 727.360   <2e-16 ***
# thetaC        1.308e+00  3.854e-03 339.390   <2e-16 ***
# tmrcaC        2.558e-02  8.061e-04  31.733   <2e-16 ***
# rhoC          6.075e-03  1.879e-02   0.323    0.747    
# thetaC:tmrcaC 1.713e+00  1.105e-01  15.505   <2e-16 ***
  
vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.036853      1.024970      1.010843      1.022137

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 9] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 9] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 9] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 9] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 9] <- anova.diversity$VarExp[4] * 100


# rep_10
rep_10.pi.200kb <- read.table("rep_10/rep_10.pi.200kb.bedgraph", header = T)
rep_10.tmrca.200kb <- read.table("rep_10/sim.tmrca.200k.map", header = T)

sim.lands.200kb <- as.data.frame(cbind(rep_10.pi.200kb$pi, sim.theta.200kb$sim, sim.rho.200kb$sim, rep_10.tmrca.200kb$tmrca))
names(sim.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.200kb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.200kb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.200kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.200kb)
plot(resid(m.div.200kb)~fitted(m.div.200kb))
dwtest(m.div.200kb)
hmctest(m.div.200kb)
hist(resid(m.div.200kb))

summary(m.div.200kb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.061e-02  2.834e-05 727.360   <2e-16 ***
# thetaC        1.308e+00  3.854e-03 339.390   <2e-16 ***
# tmrcaC        2.558e-02  8.061e-04  31.733   <2e-16 ***
# rhoC          6.075e-03  1.879e-02   0.323    0.747    
# thetaC:tmrcaC 1.713e+00  1.105e-01  15.505   <2e-16 ***

vif(m.div.200kb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.036853      1.024970      1.010843      1.022137

anova.diversity <- Anova(m.div.200kb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.200kb[1, 10] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.200kb[2, 10] <- anova.diversity$VarExp[1] * 100
r2.sim.200kb[3, 10] <- anova.diversity$VarExp[2] * 100
r2.sim.200kb[4, 10] <- anova.diversity$VarExp[3] * 100
r2.sim.200kb[5, 10] <- anova.diversity$VarExp[4] * 100


r2.sim.200kb$average <- rowMeans(r2.sim.200kb)
r2.sim.200kb <- transform(r2.sim.200kb, sd=apply(r2.sim.200kb, 1, sd, na.rm = TRUE))

# 1 Mb
r2.sim.1Mb <- as.data.frame(matrix(ncol = nreps, nrow = 5))
row.names(r2.sim.1Mb) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA")
colnames(r2.sim.1Mb) <- reps

sim.rho.1Mb <- read.table("Misc/dm.sim.rho.1e+06.txt", header = T)
sim.theta.1Mb <- read.table("Misc/dm.sim.theta.1e+06.txt", header = T)

# rep 1
rep_1.pi.1Mb <- read.table("rep_1/rep_1.pi.1Mb.bedgraph", header = T)
rep_1.tmrca.1Mb <- read.table("rep_1/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_1.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_1.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.067e-02  4.132e-05 500.181  < 2e-16 ***
# thetaC         1.313e+00  9.626e-03 136.418  < 2e-16 ***
# tmrcaC         2.602e-02  2.881e-03   9.031 2.41e-09 ***
# rhoC          -2.472e-03  6.889e-02  -0.036   0.9717    
# thetaC:tmrcaC  1.182e+00  6.258e-01   1.889   0.0706 . 

vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.248566      1.036015      1.121845      1.370034 

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 1] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 1] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 1] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 1] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 1] <- anova.diversity$VarExp[4] * 100


# rep_2
rep_2.pi.1Mb <- read.table("rep_2/rep_2.pi.1Mb.bedgraph", header = T)
rep_2.tmrca.1Mb <- read.table("rep_2/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_2.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_2.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.052e-02  4.606e-05 445.437  < 2e-16 ***
# thetaC         1.285e+00  1.013e-02 126.882  < 2e-16 ***
# tmrcaC         2.480e-02  2.820e-03   8.796 3.99e-09 ***
# rhoC          -4.630e-02  7.505e-02  -0.617 0.542900    
# thetaC:tmrcaC  2.351e+00  5.373e-01   4.376 0.000188 ***
  
vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.136589      1.161898      1.094711      1.176801 

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 2] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 2] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 2] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 2] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 2] <- anova.diversity$VarExp[4] * 100

# rep_3
rep_3.pi.1Mb <- read.table("rep_3/rep_3.pi.1Mb.bedgraph", header = T)
rep_3.tmrca.1Mb <- read.table("rep_3/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_3.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_3.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.070e-02  3.559e-05 581.752  < 2e-16 ***
# thetaC        1.326e+00  7.629e-03 173.822  < 2e-16 ***
# tmrcaC        3.071e-02  2.054e-03  14.949 5.67e-14 ***
# rhoC          6.225e-02  5.347e-02   1.164    0.255    
# thetaC:tmrcaC 2.593e+00  5.197e-01   4.990 3.83e-05 ***
  
vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC  
# 1.266893      1.266548      1.091885      1.131259 

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 3] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 3] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 3] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 3] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 3] <- anova.diversity$VarExp[4] * 100

# rep_4
rep_4.pi.1Mb <- read.table("rep_4/rep_4.pi.1Mb.bedgraph", header = T)
rep_4.tmrca.1Mb <- read.table("rep_4/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_4.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_4.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# (Intercept)   2.049e-02  3.521e-05 582.000  < 2e-16 ***
# thetaC        1.290e+00  7.811e-03 165.192  < 2e-16 ***
# tmrcaC        2.364e-02  1.989e-03  11.884 8.84e-12 ***
# rhoC          1.009e-02  5.807e-02   0.174    0.863    
# thetaC:tmrcaC 2.028e+00  3.007e-01   6.744 4.55e-07 ***

vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.146913      1.170471      1.112350      1.295650

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 4] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 4] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 4] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 4] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 4] <- anova.diversity$VarExp[4] * 100

# rep_5
rep_5.pi.1Mb <- read.table("rep_5/rep_5.pi.1Mb.bedgraph", header = T)
rep_5.tmrca.1Mb <- read.table("rep_5/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_5.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_5.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# (Intercept)    2.060e-02  3.534e-05 582.791  < 2e-16 ***
# thetaC         1.306e+00  7.849e-03 166.389  < 2e-16 ***
# tmrcaC         2.514e-02  2.157e-03  11.655 1.34e-11 ***
# rhoC          -3.495e-02  5.720e-02  -0.611   0.5467    
# thetaC:tmrcaC  1.135e+00  5.053e-01   2.246   0.0338 * 

vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.137681      1.005247      1.060240      1.119618 

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 5] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 5] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 5] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 5] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 5] <- anova.diversity$VarExp[4] * 100

# rep_6
rep_6.pi.1Mb <- read.table("rep_6/rep_6.pi.1Mb.bedgraph", header = T)
rep_6.tmrca.1Mb <- read.table("rep_6/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_6.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_6.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.055e-02  3.357e-05 612.052  < 2e-16 ***
# thetaC        1.305e+00  7.159e-03 182.260  < 2e-16 ***
# tmrcaC        2.808e-02  2.754e-03  10.193 2.18e-10 ***
# rhoC          6.955e-02  5.357e-02   1.298  0.20607    
# thetaC:tmrcaC 2.260e+00  6.739e-01   3.353  0.00255 ** 

vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.116973      1.106634      1.097630      1.052297 

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 6] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 6] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 6] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 6] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 6] <- anova.diversity$VarExp[4] * 100

# rep_7
rep_7.pi.1Mb <- read.table("rep_7/rep_7.pi.1Mb.bedgraph", header = T)
rep_7.tmrca.1Mb <- read.table("rep_7/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_7.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_7.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2.060e-02  3.433e-05 599.990  < 2e-16 ***
# thetaC         1.322e+00  7.722e-03 171.218  < 2e-16 ***
# tmrcaC         2.856e-02  1.841e-03  15.507 2.47e-14 ***
# rhoC          -5.037e-02  5.447e-02  -0.925 0.363900    
# thetaC:tmrcaC  1.507e+00  3.301e-01   4.565 0.000115 ***

vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.220441      1.073803      1.065300      1.143527

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 7] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 7] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 7] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 7] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 7] <- anova.diversity$VarExp[4] * 100


# rep_8
rep_8.pi.1Mb <- read.table("rep_8/rep_8.pi.1Mb.bedgraph", header = T)
rep_8.tmrca.1Mb <- read.table("rep_8/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_8.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_8.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.051e-02  3.685e-05 556.403  < 2e-16 ***
# thetaC        1.294e+00  7.915e-03 163.527  < 2e-16 ***
# tmrcaC        2.659e-02  2.276e-03  11.682 1.27e-11 ***
# rhoC          2.309e-02  5.876e-02   0.393  0.69767    
# thetaC:tmrcaC 2.000e+00  5.608e-01   3.567  0.00149 ** 

vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.085610      1.022656      1.049827      1.032819 

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 8] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 8] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 8] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 8] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 8] <- anova.diversity$VarExp[4] * 100

# rep_9
rep_9.pi.1Mb <- read.table("rep_9/rep_9.pi.1Mb.bedgraph", header = T)
rep_9.tmrca.1Mb <- read.table("rep_9/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_9.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_9.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.064e-02  5.033e-05 410.041  < 2e-16 ***
# thetaC        1.316e+00  1.246e-02 105.593  < 2e-16 ***
# tmrcaC        2.428e-02  2.655e-03   9.145 1.89e-09 ***
# rhoC          3.532e-02  8.105e-02   0.436   0.6667    
# thetaC:tmrcaC 1.211e+00  5.581e-01   2.171   0.0396 *  

vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.549740      1.174437      1.150241      1.464009 

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 9] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 9] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 9] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 9] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 9] <- anova.diversity$VarExp[4] * 100

# rep_10
rep_10.pi.1Mb <- read.table("rep_10/rep_10.pi.1Mb.bedgraph", header = T)
rep_10.tmrca.1Mb <- read.table("rep_10/sim.tmrca.1M.map", header = T)

sim.lands.1Mb <- as.data.frame(cbind(rep_10.pi.1Mb$pi, sim.theta.1Mb$sim, sim.rho.1Mb$sim, rep_10.tmrca.1Mb$tmrca))
names(sim.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")

cor.test(~theta+tmrca, data = sim.lands.1Mb, method = "spearman") 
cor.test(~theta+rho, data = sim.lands.1Mb, method = "spearman") 
cor.test(~rho+tmrca, data = sim.lands.1Mb, method = "spearman") 

cor.test(~diversity+theta, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+tmrca, data = sim.lands.200kb, method = "spearman") 
cor.test(~diversity+rho, data = sim.lands.200kb, method = "spearman") 

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

m.div.1Mb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = sim.lands.1Mb)
plot(resid(m.div.1Mb)~fitted(m.div.1Mb))
dwtest(m.div.1Mb)
hmctest(m.div.1Mb)
hist(resid(m.div.1Mb))

summary(m.div.1Mb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.064e-02  5.033e-05 410.041  < 2e-16 ***
# thetaC        1.316e+00  1.246e-02 105.593  < 2e-16 ***
# tmrcaC        2.428e-02  2.655e-03   9.145 1.89e-09 ***
# rhoC          3.532e-02  8.105e-02   0.436   0.6667    
# thetaC:tmrcaC 1.211e+00  5.581e-01   2.171   0.0396 *  

vif(m.div.1Mb)
# thetaC        tmrcaC          rhoC thetaC:tmrcaC 
# 1.549740      1.174437      1.150241      1.464009 

anova.diversity <- Anova(m.div.1Mb)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

r2.sim.1Mb[1, 10] <- (anova.diversity$VarExp[1] + anova.diversity$VarExp[2] + anova.diversity$VarExp[3] + anova.diversity$VarExp[4]) * 100
r2.sim.1Mb[2, 10] <- anova.diversity$VarExp[1] * 100
r2.sim.1Mb[3, 10] <- anova.diversity$VarExp[2] * 100
r2.sim.1Mb[4, 10] <- anova.diversity$VarExp[3] * 100
r2.sim.1Mb[5, 10] <- anova.diversity$VarExp[4] * 100


r2.sim.1Mb$average <- rowMeans(r2.sim.1Mb)
r2.sim.1Mb <- transform(r2.sim.1Mb, sd=apply(r2.sim.1Mb, 1, sd, na.rm = TRUE))



#######################################
#
# R2 plot
#
#######################################

# true landscapes
r2.sim.avg <- rbind.data.frame(r2.sim.50kb$average, r2.sim.200kb$average, r2.sim.1Mb$average, make.row.names = F)
colnames(r2.sim.avg) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA")
r2.sim.avg$bin.size <- c(50, 200, 1000)

# inferred landscapes
r2.inf.avg <- rbind.data.frame(r2.inf.50kb$average, r2.inf.200kb$average, r2.inf.1Mb$average, make.row.names = F)
colnames(r2.inf.avg) <- c("Total", "Theta", "Rho", "TMRCA", "Theta:TMRCA")
r2.inf.avg$bin.size <- c(50, 200, 1000)
  
r2.dm.tab.2 <- as.data.frame(cbind(apply(r2.dm.tab, 2, as.numeric)))

r2.tab.comb <- rbind.data.frame(r2.tab.2, r2.sim.avg, r2.inf.avg)
r2.tab.comb$type <- c(rep("real", 3), rep("sim", 3), rep("inf", 3))

molten.r2 <- melt(r2.tab.comb, id.vars = c("bin.size", "type"))
r2.plot <- ggplot(data = molten.r2, aes(x = bin.size, y = value, colour = variable, fill = type))
r2.plot <- r2.plot + geom_line(data = molten.r2)
r2.plot <- r2.plot + geom_point(aes(shape = type, colour = variable), size = 4)
r2.plot <- r2.plot + scale_x_continuous(breaks = c(50, 200, 1000)) 
r2.plot <- r2.plot + scale_y_continuous(breaks = pretty_breaks())
r2.plot <- r2.plot + labs(title = NULL, x = "Bin Size (kb)", y = "Variance Explained (%)") + theme.blank
r2.plot <- r2.plot + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16))

leg <- get_legend(r2.plot + theme(legend.position="bottom"))

fig4 <- plot_grid(r2.plot + no.legend, leg, nrow = 2, labels = NULL, rel_heights = c(1, 0.1), scale = c(1, 0.2))

cowplot::save_plot("../submission/Figure4.pdf", plot = fig4, base_width = 10, base_height = 6, device = "pdf")

