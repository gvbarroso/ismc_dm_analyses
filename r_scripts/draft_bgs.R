

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

setwd("~/Data/iSMC/theta_paper/BGS_iSMC/30x5x5/")

# for plotting exons later on:
dm_2L_exome <- read.table("~/Data/iSMC/theta_paper/slim_sims/dm_tbl.txt", header = T) %>% 
  filter(seqnames == "1") %>%
  filter(!is.na(ensembl_exon_id)) %>%
  dplyr::select(start, end)

names(dm_2L_exome) <- c("exon_start", "exon_end")

bgs.TMRCA.50kb <- read.table("dm_2L_bgs_rep_1.TMRCA.50kb.bedgraph", header = T)
bgs.TMRCA.200kb <- read.table("dm_2L_bgs_rep_1.TMRCA.200kb.bedgraph", header = T)
bgs.TMRCA.1Mb <- read.table("dm_2L_bgs_rep_1.TMRCA.1Mb.bedgraph", header = T)

bgs.TMRCA.50kb$mean <- apply(bgs.TMRCA.50kb[,(4:ncol(bgs.TMRCA.50kb))], 1, mean)
bgs.TMRCA.200kb$mean <- apply(bgs.TMRCA.200kb[,(4:ncol(bgs.TMRCA.200kb))], 1, mean)
bgs.TMRCA.1Mb$mean <- apply(bgs.TMRCA.1Mb[,(4:ncol(bgs.TMRCA.1Mb))], 1, mean)

bgs.rho.50kb <- read.table("dm_2L_bgs_rep_1.rho.50kb.bedgraph", header = T)
bgs.rho.200kb <- read.table("dm_2L_bgs_rep_1.rho.200kb.bedgraph", header = T)
bgs.rho.1Mb <- read.table("dm_2L_bgs_rep_1.rho.1Mb.bedgraph", header = T)

bgs.theta.50kb <- read.table("dm_2L_bgs_rep_1.theta.50kb.bedgraph", header = T)
bgs.theta.200kb <- read.table("dm_2L_bgs_rep_1.theta.200kb.bedgraph", header = T)
bgs.theta.1Mb <- read.table("dm_2L_bgs_rep_1.theta.1Mb.bedgraph", header = T)

bgs.pi.50kb <- read.table("dm_2L_bgs_rep_1.diversity.50kb.bedgraph", header = T)
bgs.pi.200kb <- read.table("dm_2L_bgs_rep_1.diversity.200kb.bedgraph", header = T)
bgs.pi.1Mb <- read.table("dm_2L_bgs_rep_1.diversity.1Mb.bedgraph", header = T)

bgs.pi.50kb$mean <- apply(bgs.pi.50kb[,(4:ncol(bgs.pi.50kb))], 1, mean)
bgs.pi.200kb$mean <- apply(bgs.pi.200kb[,(4:ncol(bgs.pi.200kb))], 1, mean)
bgs.pi.1Mb$mean <- apply(bgs.pi.1Mb[,(4:ncol(bgs.pi.1Mb))], 1, mean)

# 50kb
bgs.lands.50kb <- as.data.frame(cbind(bgs.pi.50kb$mean, bgs.theta.50kb$sample_mean, bgs.rho.50kb$sample_mean, bgs.TMRCA.50kb$mean))
names(bgs.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")

# centering
bgs.lands.50kb$thetaC <- bgs.lands.50kb$theta - mean(bgs.lands.50kb$theta)
bgs.lands.50kb$tmrcaC <- bgs.lands.50kb$tmrca - mean(bgs.lands.50kb$tmrca)
bgs.lands.50kb$rhoC <- bgs.lands.50kb$rho - mean(bgs.lands.50kb$rho)

bgs.lands.50kb$bin <- 1:nrow(bgs.lands.50kb)

m.bgs.50kb <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC, data = bgs.lands.50kb)
m.bgs.50kb.2 <- lm(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC + rhoC:tmrcaC, data = bgs.lands.50kb)
m.bgs.50kb.3 <- lm(diversity ~ (thetaC + rhoC + tmrcaC) ^ 2, data = bgs.lands.50kb)

AIC(m.bgs.50kb, m.bgs.50kb.2, m.bgs.50kb.3)

plot(resid(m.bgs.50kb.3)~fitted(m.bgs.50kb.3))
dwtest(m.bgs.50kb.3)
hmctest(m.bgs.50kb.3)
hist(resid(m.bgs.50kb.3))

summary(m.bgs.50kb.3)

anova.diversity <- Anova(m.bgs.50kb.3)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

#Anova Table (Type II tests)

#Response: diversity
#Sum Sq  Df    F value  Pr(>F)  VarExp
#thetaC        2.6580e-07   1  1017.4294 0.00000 0.04959
#rhoC          4.0000e-10   1     1.6387 0.20114 0.00008
#tmrcaC        4.8618e-06   1 18611.2157 0.00000 0.90711
#thetaC:rhoC   2.0000e-10   1     0.6039 0.43747 0.00003
#thetaC:tmrcaC 1.1060e-07   1   423.5635 0.00000 0.02064
#rhoC:tmrcaC   2.0000e-10   1     0.6500 0.42054 0.00003
#Residuals     1.2070e-07 462                    0.02252

g.bgs.50kb.1 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = bgs.lands.50kb, weights = varPower(0, ~tmrcaC), cor = corAR1(0, ~bin), method = "ML")

g.bgs.50kb.2 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = bgs.lands.50kb, weights = varPower(0, ~theta), cor = corAR1(0, ~bin), method = "ML")

g.bgs.50kb.3 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = bgs.lands.50kb, weights = varPower(0, ~theta), method = "ML")

g.bgs.50kb.4 <- gls(diversity ~ thetaC + rhoC + tmrcaC + thetaC:tmrcaC,
                    data = bgs.lands.50kb, weights = varPower(0, ~tmrcaC), method = "ML")

AIC(g.bgs.50kb.1, g.bgs.50kb.2, g.bgs.50kb.3, g.bgs.50kb.4)

summary(g.bgs.50kb.1)
#Coefficients:
#  Value  Std.Error   t-value p-value
#(Intercept)    0.000117 0.00000070 167.74820  0.0000
#thetaC         5.174353 0.20967174  24.67835  0.0000
#rhoC          -0.014365 0.01626844  -0.88300  0.3777
#tmrcaC         0.000145 0.00000117 124.15330  0.0000
#thetaC:tmrcaC  4.080465 0.21313134  19.14531  0.0000

vif(g.bgs.50kb.3)
#thetaC          rhoC        tmrcaC thetaC:tmrcaC 
#1.072982      1.072728      1.354150      1.374055 

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


