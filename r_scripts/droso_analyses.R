# Created: 03/07/2019
# Last modified: 11/08/2020
# Author: Gustavo Barroso

library(ppcor)
library(MASS)
library(reshape2)
library(ggplot2)
library(ppls)
library(scales)
library(cowplot)
library(lmtest)
library(relaimpo)
library(corrplot)
library(nlme)
library(car)
library(plyr)
library(GenomicRanges)

setwd("~")
setwd("Data/iSMC/theta_paper/real_data/droso.chr2L.phased/")

# R_2 table for plotting at the end
r2.tab <- as.data.frame(matrix(ncol = 5, nrow = 3))
colnames(r2.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Bin_size(kb)")

theme.blank <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

######################################
#
# Drosophila-like neutral simulations T
# TODO
#
########################################

# building linear models

nreps <- 10

# 50 kb
for(in in 1:nreps)
{
  p <- paste("bedgraph/rs.pair_", i, ".", sep = "")
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
  anova.diversity <- Anova(m.diversity.bc)
  apiss <- anova.diversity$"Sum Sq"
  anova.diversity$VarExp <- apiss / sum(apiss)
  
  anova.diversity
}


######################################
#
# 30x5x5 --- 50kb
#
########################################

# recombination landscapes
rho.dm.50kb <- read.table("dm_30x5x5.rho.50kb.bedgraph", header = T)

# diversity
diversity.dm.50kb <- read.table("dm_30x5x5.diversity.50kb.bedgraph", header = T)
diversity.dm.50kb$avg <- apply(diversity.dm.50kb[4:ncol(diversity.dm.50kb)], 1, mean)

# mutation landscapes
theta.dm.50kb <- read.table("dm_30x5x5.theta.50kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.50kb <- read.table("dm_30x5x5.TMRCA.50kb.bedgraph", header = T)
tmrca.dm.50kb$sample_mean <- apply(tmrca.dm.50kb[4:ncol(tmrca.dm.50kb)], 1, mean)

# missing data
missing.prop.50kb <- read.table("dm_30x5x5.missing.prop.50kb.bedgraph", header = T)
intersect.50kb <- apply(missing.prop.50kb[4:ncol(missing.prop.50kb)], 1, function(x) any(x > 0.25)) 

dm.lands.50kb <- as.data.frame(cbind(diversity.dm.50kb$avg, theta.dm.50kb$sample_mean, rho.dm.50kb$sample_mean, tmrca.dm.50kb$sample_mean))
names(dm.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
dm.lands.50kb$bin <- 1:nrow(dm.lands.50kb)
# filters based on missing data ( > 50% per window)
dm.lands.50kb <- dm.lands.50kb[which(intersect.50kb == F),]

# multi-collinearity, unlike neutral simulations
cor.mat <- cor(dm.lands.50kb[,1:4], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# we try to disentangle using partial correlations
pcor.mat <- cor.mat

pcor.test(x = dm.lands.50kb$rho, y = dm.lands.50kb$tmrca, z = dm.lands.50kb$theta, method = "spearman")
pcor.mat [4, 3] <- 0.40
pcor.test(x = dm.lands.50kb$rho, y = dm.lands.50kb$diversity, z = dm.lands.50kb$tmrca, method = "spearman")
pcor.test(x = dm.lands.50kb$rho, y = dm.lands.50kb$diversity, z = dm.lands.50kb$theta, method = "spearman")
pcor.test(x = dm.lands.50kb$rho, y = dm.lands.50kb$diversity, z = cbind(dm.lands.50kb$tmrca, dm.lands.50kb$theta), method = "spearman")
pcor.mat[3, 1] <- 0
pcor.test(x = dm.lands.50kb$rho, y = dm.lands.50kb$theta, z = dm.lands.50kb$tmrca, method = "spearman")
pcor.mat[3, 2] <- 0

pcor.test(x = dm.lands.50kb$theta, y = dm.lands.50kb$diversity, z = dm.lands.50kb$tmrca, method = "spearman")
pcor.test(x = dm.lands.50kb$theta, y = dm.lands.50kb$diversity, z = dm.lands.50kb$rho, method = "spearman")
pcor.test(x = dm.lands.50kb$theta, y = dm.lands.50kb$diversity, z = cbind(dm.lands.50kb$rho, dm.lands.50kb$tmrca), method = "spearman")
pcor.mat[2, 1] <- 0.97

pcor.test(x = dm.lands.50kb$theta, y = dm.lands.50kb$tmrca, z = dm.lands.50kb$rho, method = "spearman")
pcor.mat[4, 2] <- 0.48

pcor.test(x = dm.lands.50kb$tmrca, y = dm.lands.50kb$diversity, z = dm.lands.50kb$theta, method = "spearman")
pcor.test(x = dm.lands.50kb$tmrca, y = dm.lands.50kb$diversity, z = cbind(dm.lands.50kb$theta, dm.lands.50kb$rho), method = "spearman")
pcor.mat[4, 1] <- 0.72

pdf("dm.pcor.50kb.pdf")
corrplot(pcor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()


# CV (added to plots with annotate)
sd(dm.lands.50kb$diversity) / mean(dm.lands.50kb$diversity) # 0.27
sd(dm.lands.50kb$theta) / mean(dm.lands.50kb$theta) # 0.24
sd(dm.lands.50kb$rho) / mean(dm.lands.50kb$rho) # 0.22
sd(dm.lands.50kb$tmrca) / mean(dm.lands.50kb$tmrca) # 0.07

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(dm.lands.50kb[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 20, y = value)) 
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
diversity.map <- diversity.map + annotate("text", x = 4100, y = 0.005, label = "CV = 27%")

molten.rho <- melt(dm.lands.50kb[c(3,5)], id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 20, y = value)) 
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
rho.map <- rho.map + annotate("text", x = 4100, y = 0.06, label = "CV = 22%")

molten.theta <- melt(dm.lands.50kb[c(2,5)], id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 20, y = value)) 
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
theta.map <- theta.map + annotate("text", x = 4100, y = 0.0145, label = "CV = 24%")

molten.tmrca <- melt(dm.lands.50kb[c(4,5)], id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 20, y = value)) 
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
tmrca.map <- tmrca.map + annotate("text", x = 4100, y = 0.7, label = "CV = 7%")

p1 <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1)
p1

cowplot::save_plot("dm.maps.50kb.pdf", plot = p, device = "pdf", dpi = 500, base_width = 10, base_height = 8)




# Linear model with all variables
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.50kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) #***
hist(resid(m.diversity.bc), nclass = 30) # not too bad?
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # ***
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # **

summary(m.diversity.bc)
# (Intercept) -1.4901599  0.0008151 -1828.172   <2e-16 ***
# theta        4.4862107  0.0267419   167.760   <2e-16 ***
# rho          0.0049350  0.0079687     0.619    0.536    
# tmrca        0.0498502  0.0010696    46.605   <2e-16 ***
# Multiple R-squared: 0.9923,	Adjusted R-squared:  0.9923 

vif(m.diversity.bc)
# theta      rho    tmrca 
# 1.335152 1.382633 1.755700 

# type 2
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# theta     0.032771   1 28143.2513 0.00000 0.91586
# rho       0.000000   1     0.3835 0.53606 0.00001
# tmrca     0.002529   1  2172.0638 0.00000 0.07069
# Residuals 0.000481 413                    0.01344

r2.tab[1,] <- c(91.6 + 7.1, 91.6, 0.0, 7.1, 50)

# Because of auto-correlation we compute p-values for the variables using a GLS
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = dm.lands.50kb, corr = corAR1(0, ~bin))

summary(g.diversity)
# Parameter estimate(s):
#  Phi 
# 0.3903864 

# Coefficients:
# (Intercept) -0.0109222 0.000247337 -44.15935   0.000
# theta        0.9656258 0.010592500  91.16127   0.000
# rho          0.0035266 0.002614014   1.34910   0.178
# tmrca        0.0118165 0.000322324  36.66026   0.000



# Linear model without TMRCA  -> rho becomes significant
m.diversity <- lm(diversity ~ rho + theta, data = dm.lands.50kb)
plot(m.diversity, which = 2)
summary(m.diversity)
# (Intercept) -0.0023976  0.0001694  -14.15   <2e-16 ***
# rho          0.0435202  0.0038502   11.30   <2e-16 ***
# theta        1.0974777  0.0131487   83.47   <2e-16 ***
# Multiple R-squared:  0.9499,	Adjusted R-squared:  0.9496 



########################################
#
# 30x5x5 --- 200kb
#
########################################

# recombination landscapes
rho.dm.200kb <- read.table("dm_30x5x5.rho.200kb.bedgraph", header = T)

# diversity
diversity.dm.200kb <- read.table("dm_30x5x5.diversity.200kb.bedgraph", header = T)
diversity.dm.200kb$avg <- apply(diversity.dm.200kb[4:ncol(diversity.dm.200kb)], 1, mean)

# mutation landscapes
theta.dm.200kb <- read.table("dm_30x5x5.theta.200kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.200kb <- read.table("dm_30x5x5.TMRCA.200kb.bedgraph", header = T)
tmrca.dm.200kb$sample_mean <- apply(tmrca.dm.200kb[4:ncol(tmrca.dm.200kb)], 1, mean)

# missing data
missing.prop.200kb <- read.table("dm_30x5x5.missing.prop.200kb.bedgraph", header = T)
intersect.200kb <- apply(missing.prop.200kb[4:ncol(missing.prop.200kb)], 1, function(x) any(x > 0.25)) 

dm.lands.200kb <- as.data.frame(cbind(diversity.dm.200kb$avg, theta.dm.200kb$sample_mean, rho.dm.200kb$sample_mean, tmrca.dm.200kb$sample_mean))
names(dm.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
dm.lands.200kb$bin <- 1:nrow(dm.lands.200kb)
# filters based on missing data ( > 50% per window)
dm.lands.200kb <- dm.lands.200kb[which(intersect.200kb == F),]

# multi-collinearity, unlike neutral simulations
cor.mat <- cor(dm.lands.200kb[,1:4], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# we try to disentangle using partial correlations
pcor.mat <- cor.mat

pcor.test(x = dm.lands.200kb$rho, y = dm.lands.200kb$tmrca, z = dm.lands.200kb$theta, method = "spearman")
pcor.mat[4, 3] <- 0.49
pcor.test(x = dm.lands.200kb$rho, y = dm.lands.200kb$diversity, z = dm.lands.200kb$tmrca, method = "spearman")
pcor.mat[3, 1] <- 0
pcor.test(x = dm.lands.200kb$rho, y = dm.lands.200kb$theta, z = dm.lands.200kb$tmrca, method = "spearman")
pcor.mat[3, 2] <- 0

pcor.test(x = dm.lands.200kb$theta, y = dm.lands.200kb$diversity, z = dm.lands.200kb$tmrca, method = "spearman")
pcor.test(x = dm.lands.200kb$theta, y = dm.lands.200kb$diversity, z = dm.lands.200kb$rho, method = "spearman")
pcor.test(x = dm.lands.200kb$theta, y = dm.lands.200kb$diversity, z = cbind(dm.lands.200kb$rho, dm.lands.200kb$tmrca), method = "spearman")
pcor.mat[2, 1] <- 0.98
pcor.test(x = dm.lands.200kb$theta, y = dm.lands.200kb$tmrca, z = dm.lands.200kb$rho, method = "spearman")
pcor.mat[4, 2] <- 0.52

pcor.test(x = dm.lands.200kb$tmrca, y = dm.lands.200kb$diversity, z = dm.lands.200kb$rho, method = "spearman")
pcor.test(x = dm.lands.200kb$tmrca, y = dm.lands.200kb$diversity, z = dm.lands.200kb$theta, method = "spearman")
pcor.test(x = dm.lands.200kb$tmrca, y = dm.lands.200kb$diversity, z = cbind(dm.lands.200kb$theta, dm.lands.200kb$rho), method = "spearman")
pcor.mat[4, 1] <- 0.65

pdf("dm.pcor.200kb.pdf")
corrplot(pcor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()

# CV (added to plots with annotate)
sd(dm.lands.200kb$diversity) / mean(dm.lands.200kb$diversity) # 0.24
sd(dm.lands.200kb$theta) / mean(dm.lands.200kb$theta) # 0.21
sd(dm.lands.200kb$rho) / mean(dm.lands.200kb$rho) # 0.17
sd(dm.lands.200kb$tmrca) / mean(dm.lands.200kb$tmrca) # 0.05

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(dm.lands.200kb[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 200, y = value)) 
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
diversity.map <- diversity.map + annotate("text", x = 11000, y = 0.005, label = "CV = 24%")

molten.rho <- melt(dm.lands.200kb[c(3,5)], id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 200, y = value)) 
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 50), axis.title.y = element_text(size = 20)) + theme.blank
rho.map <- rho.map + annotate("text", x = 11000, y = 0.0475, label = "CV = 16%")

molten.theta <- melt(dm.lands.200kb[c(2,5)], id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 200, y = value)) 
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
theta.map <- theta.map + annotate("text", x = 11000, y = 0.006, label = "CV = 21%")

molten.tmrca <- melt(dm.lands.200kb[c(4,5)], id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 200, y = value)) 
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
tmrca.map <- tmrca.map + annotate("text", x = 11000, y = 0.75, label = "CV = 5%")

p2 <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1)
p2


cowplot::save_plot("dm.maps.200kb.pdf", plot = p, device = "pdf", dpi = 500, base_width = 10, base_height = 8)





# Linear model with all variables
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.200kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # NS
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # *
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# (Intercept) -1.346506   0.001119 -1203.475   <2e-16 ***
# theta        3.229108   0.028949   111.545   <2e-16 ***
# rho          0.014543   0.009857     1.475    0.143    
# tmrca        0.031941   0.001494    21.381   <2e-16 ***
# Multiple R-squared:  0.9961,	Adjusted R-squared:  0.996

vif(m.diversity.bc)
# theta      rho    tmrca 
# 1.559019 1.484401 2.151090 

# type 2 anova
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# theta     0.00297305   1 12442.3642 0.00000 0.95691
# rho       0.00000052   1     2.1767 0.14322 0.00017
# tmrca     0.00010923   1   457.1430 0.00000 0.03516
# Residuals 0.00002413 101                    0.00777

r2.tab[2,] <- c(95.7 + 3.5, 95.7, 0, 3.5, 200)

# Because of auto-correlation we compute p-values for the variables using a GLS
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = dm.lands.200kb, corr = corAR1(0, ~bin))

summary(g.diversity)
# Parameter estimate(s):
#  Phi 
# 0.2781892

# Coefficients:
#  Value   Std.Error   t-value p-value
# (Intercept) -0.0080258 0.000461083 -17.40640  0.0000
# theta        0.9956466 0.013381379  74.40538  0.0000
# rho          0.0125460 0.004487122   2.79601  0.0062
# tmrca        0.0081182 0.000621129  13.07013  0.0000




########################################
#
# 30x5x5 --- 1Mb
#
########################################

# recombination landscapes
rho.dm.1Mb <- read.table("dm_30x5x5.rho.1Mb.bedgraph", header = T)

# diversity
diversity.dm.1Mb <- read.table("dm_30x5x5.diversity.1Mb.bedgraph", header = T)
diversity.dm.1Mb$avg <- apply(diversity.dm.1Mb[4:ncol(diversity.dm.1Mb)], 1, mean)

# mutation landscapes
theta.dm.1Mb <- read.table("dm_30x5x5.theta.1Mb.bedgraph", header = T)

# TMRCA landscapes
tmrca.dm.1Mb <- read.table("dm_30x5x5.TMRCA.1Mb.bedgraph", header = T)
tmrca.dm.1Mb$sample_mean <- apply(tmrca.dm.1Mb[4:ncol(tmrca.dm.1Mb)], 1, mean)

# missing data
missing.prop.1Mb <- read.table("dm_30x5x5.missing.prop.1Mb.bedgraph", header = T)
intersect.1Mb <- apply(missing.prop.1Mb[4:ncol(missing.prop.1Mb)], 1, function(x) any(x > 0.25)) 

dm.lands.1Mb <- as.data.frame(cbind(diversity.dm.1Mb$avg, theta.dm.1Mb$sample_mean, rho.dm.1Mb$sample_mean, tmrca.dm.1Mb$sample_mean))
names(dm.lands.1Mb) <- c("diversity", "theta", "rho", "tmrca")
dm.lands.1Mb$bin <- 1:nrow(dm.lands.1Mb)
# filters based on missing data ( > 50% per window)
dm.lands.1Mb <- dm.lands.1Mb[which(intersect.1Mb == F),]

# multi-collinearity, unlike neutral simulations
cor.mat <- cor(dm.lands.1Mb[,1:4], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# we try to disentangle using partial correlations
pcor.mat <- cor.mat 

pcor.test(x = dm.lands.1Mb$rho, y = dm.lands.1Mb$tmrca, z = dm.lands.1Mb$theta, method = "spearman")
pcor.mat[4, 3] <- 0.52
pcor.test(x = dm.lands.1Mb$rho, y = dm.lands.1Mb$diversity, z = dm.lands.1Mb$tmrca, method = "spearman")
pcor.mat[3, 1] <- 0
pcor.test(x = dm.lands.1Mb$rho, y = dm.lands.1Mb$theta, z = dm.lands.1Mb$tmrca, method = "spearman")
pcor.mat[3, 2] <- 0

pcor.test(x = dm.lands.1Mb$theta, y = dm.lands.1Mb$diversity, z = dm.lands.1Mb$tmrca, method = "spearman")
pcor.test(x = dm.lands.1Mb$theta, y = dm.lands.1Mb$diversity, z = dm.lands.1Mb$rho, method = "spearman")
pcor.test(x = dm.lands.1Mb$theta, y = dm.lands.1Mb$diversity, z = cbind(dm.lands.1Mb$rho, dm.lands.1Mb$tmrca), method = "spearman")
pcor.mat[2, 1] <- 0.96
pcor.test(x = dm.lands.1Mb$theta, y = dm.lands.1Mb$tmrca, z = dm.lands.1Mb$rho, method = "spearman")
pcor.mat[4, 2] <- 0.77

pcor.test(x = dm.lands.1Mb$tmrca, y = dm.lands.1Mb$diversity, z = dm.lands.1Mb$rho, method = "spearman")
pcor.test(x = dm.lands.1Mb$tmrca, y = dm.lands.1Mb$diversity, z = dm.lands.1Mb$theta, method = "spearman")
pcor.test(x = dm.lands.1Mb$tmrca, y = dm.lands.1Mb$diversity, z = cbind(dm.lands.1Mb$rho, dm.lands.1Mb$theta), method = "spearman")
pcor.mat[4, 1] <- 0

pdf("dm.pcor.1Mb.pdf")
corrplot(pcor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()

# CV (added to plots with annotate)
sd(dm.lands.1Mb$diversity) / mean(dm.lands.1Mb$diversity) # 0.16
sd(dm.lands.1Mb$theta) / mean(dm.lands.1Mb$theta) # 0.15
sd(dm.lands.1Mb$rho) / mean(dm.lands.1Mb$rho) # 0.12
sd(dm.lands.1Mb$tmrca) / mean(dm.lands.1Mb$tmrca) # 0.03

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(dm.lands.1Mb[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 1000, y = value)) 
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
diversity.map <- diversity.map + annotate("text", x = 11000, y = 0.008, label = "CV = 16%")

molten.rho <- melt(dm.lands.1Mb[c(3,5)], id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 1000, y = value)) 
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
rho.map <- rho.map + annotate("text", x = 11000, y = 0.0425, label = "CV = 12%")

molten.theta <- melt(dm.lands.1Mb[c(2,5)], id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 1000, y = value)) 
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
theta.map <- theta.map + annotate("text", x = 11000, y = 0.008, label = "CV = 15%")

molten.tmrca <- melt(dm.lands.1Mb[c(4,5)], id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 1000, y = value)) 
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) + theme.blank
tmrca.map <- tmrca.map + annotate("text", x = 11000, y = 0.9, label = "CV = 3%")

p3 <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1)
p3

cowplot::save_plot("dm.maps.1Mb.pdf", plot = p, device = "pdf", dpi = 500, base_width = 10, base_height = 8)




# Linear model with all variables
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.1Mb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # NS
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # NS
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# (Intercept) -2.150321   0.009875 -217.764  < 2e-16 ***
# theta       12.279705   0.220323   55.735  < 2e-16 ***
# rho          0.087535   0.056389    1.552     0.14    
# tmrca        0.125593   0.012563    9.997 2.76e-08 ***
# Multiple R-squared:  0.9982,	Adjusted R-squared:  0.9979 

# type 2 anova
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# theta     2.867e-03  1 3106.3864 0.00000 0.96330
# rho       2.220e-06  1    2.4097 0.14014 0.00075
# tmrca     9.224e-05  1   99.9418 0.00000 0.03099
# Residuals 1.477e-05 16                   0.00496

r2.tab[3,] <- c(96.3 + 3.1, 96.3, 0, 3.1, 1000)

# Because of auto-correlation we compute p-values for the variables using a GLS
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = dm.lands.1Mb, corr = corAR1(0, ~bin))

summary(g.diversity)
# Parameter estimate(s):
#  Phi 
# 0.2781892

# Coefficients:
#  Value   Std.Error  t-value p-value
# (Intercept) -0.0083513 0.001171002 -7.13176  0.0000
# theta        0.9887601 0.025653413 38.54302  0.0000
# rho          0.0170698 0.006551188  2.60560  0.0191
# tmrca        0.0083624 0.001481556  5.64431  0.0000


########################################
#
# Maps Plot [WEIRD ISSUE WITH FONT SIZE FOR P1]
#
########################################


pcomb <- plot_grid(p1, p2, p3, nrow = 1, labels = "AUTO", label_size = 18, scale = 0.9)
cowplot::save_plot("dm.maps.pdf", pcomb, base_height = 6, base_width = 15)

########################################
#
# R2 Plot
#
########################################

r2.tab.2 <- as.data.frame(cbind(apply(r2.tab, 2, as.numeric)))
names(r2.tab.2)[5] <- "bin.size"

molten.r2 <- melt(r2.tab.2, id.vars = "bin.size")
r2.plot <- ggplot(data = molten.r2, aes(x = bin.size, y = value, colour = variable))
r2.plot <- r2.plot + geom_line(data = molten.r2)
r2.plot <- r2.plot + geom_point(aes(colour = variable), size = 7, shape = 19)
r2.plot <- r2.plot + scale_x_continuous(breaks = c(50, 200, 1000)) 
r2.plot <- r2.plot + scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 100))
r2.plot <- r2.plot + labs(title = NULL, x = "Bin Size (kb)", y = "Var. Explained")
r2.plot <- r2.plot + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

ggsave("lm.r2.dm.pdf", plot = r2.plot, dpi = 500, width = 7, height = 7)


########################################
#
# Divergence --- 50 kb
#
########################################

# divergence data from D. melanogaster and D. yakuba
div <- read.table("~/Data/iSMC/theta_paper/real_data/droso.misc/Droso2L_divergence.statistics5kb.csv", header = T)
div <- div[which(div$Chr == "2L"), c(1:3, 6)] 

# to get the ranges
diversity.dm <- diversity.dm.50kb[which(intersect.50kb == F),] 
dm.50kb <- as.data.frame(cbind(diversity.dm[,(1:3)], dm.lands.50kb[,1:4]))
names(dm.50kb) <- c("Chr", "Start", "Stop", names(dm.lands.50kb)[1:4])
dm.50kb$Chr <- "2L"

# converts objects
dm.gr <- makeGRangesFromDataFrame(dm.50kb)
values(dm.gr) <- dm.50kb[,(4:7)]
div.gr <- makeGRangesFromDataFrame(div)
values(div.gr) <- DataFrame(score = div$MLModelFit.BrLen0)

hits <- findOverlaps(dm.gr, div.gr) 
ranges(div.gr)[subjectHits(hits)] = ranges(dm.gr)[queryHits(hits)]

lands.gr.df <- as.data.frame(dm.gr)
div.gr.df <- as.data.frame(div.gr)
# for some reason GRanges has problems converting start coordinates within 6 lands bins
div.gr.df <- div.gr.df[-which(((div.gr.df$width - 1) %% 50000) != 0),]

tmp <- ddply(.data = div.gr.df[-c(1,5)], .variables = "start", .fun = colMeans, .progress = "text")
div.gr.df.2 <- div.gr.df[!duplicated(div.gr.df$start),]
div.gr.df.2$score <- tmp$score

lands.gr.df.2 <- lands.gr.df[which(lands.gr.df$start %in% div.gr.df.2$start),]

lands.div.dm <- cbind.data.frame(lands.gr.df.2[,-(which(names(lands.gr.df.2) == "strand"))], div.gr.df.2$score)
names(lands.div.dm)[ncol(lands.div.dm)] <- "divergence"

cor.test(lands.div.dm$theta, lands.div.dm$divergence, method = "spearman")
# S = 1672646, p-value = 1.796e-05
# sample estimates:
#   rho 
# 0.2740126 

# NOTE cannot build a linear model due to extreme multi-collinearity (huge VIFs -- not shown)
cor.mat <- cor(lands.div.dm[,5:9], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# Moving on to partial correlations...

# NS
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$tmrca, z = cbind(lands.div.dm$theta, lands.div.dm$rho), method = "spearman")
# NS
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$rho, z = cbind(lands.div.dm$theta, lands.div.dm$tmrca), method = "spearman")

# afaik there are two hypothesis for the correlation between divergence and diversity:
# first, mutation rate variation
# second, background selection
cor.test(x = lands.div.dm$divergence, y = lands.div.dm$diversity, method = "spearman")
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$diversity, lands.div.dm$theta, method = "spearman")
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$diversity, z = cbind(lands.div.dm$theta, lands.div.dm$rho), method = "spearman")

# NS
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$tmrca, lands.div.dm$theta, method = "spearman")
# 0.17, p = 0.011
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$theta, z = cbind(lands.div.dm$rho, lands.div.dm$tmrca), method = "spearman")

# plot
theta.div.scatter <- ggplot(data = lands.div.dm, aes(x = lands.div.dm$theta, y = lands.div.dm$divergence, color = lands.div.dm$tmrca))
theta.div.scatter <- theta.div.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
#theta.div.scatter <- theta.div.scatter + geom_smooth(method = "loess", colour = 'magenta', linetype = 2, se = F)
theta.div.scatter <- theta.div.scatter + theme(plot.title = element_text(size = 20),
                                               axis.title.x = element_text(size = 20),
                                               axis.title.y = element_text(size = 20))
theta.div.scatter <- theta.div.scatter + labs(title = NULL, # "Mutation rate vs Divergence D. melanogaster : D. yakuba",
                                              color = "TMRCA",
                                              x = "Scaled Mutation Rate in D. melanogaster",
                                              y = "Divergence")
theta.div.scatter <- theta.div.scatter + annotate("text", x = 0.008, y = 0.095, label = "Spearman's rank partial correlation = 0.17, p = 0.011")

ggsave("divergence.dm.50kb.pdf", plot = theta.div.scatter, dpi = 500, width = 7, height = 7)


########################################
#
# Divergence --- 200 kb
#
########################################

# divergence data from D. melanogaster and D. yakuba
div <- read.table("~/Data/iSMC/theta_paper/real_data/droso.misc/Droso2L_divergence.statistics5kb.csv", header = T)
div <- div[which(div$Chr == "2L"), c(1:3, 6)] 

# to get the ranges
diversity.dm <- diversity.dm.200kb[which(intersect.200kb == F),] 
dm.200kb <- as.data.frame(cbind(diversity.dm[,(1:3)], dm.lands.200kb[,1:4]))
names(dm.200kb) <- c("Chr", "Start", "Stop", names(dm.lands.200kb)[1:4])
dm.200kb$Chr <- "2L"

# converts objects
dm.gr <- makeGRangesFromDataFrame(dm.200kb)
values(dm.gr) <- dm.200kb[,(4:7)]
div.gr <- makeGRangesFromDataFrame(div)
values(div.gr) <- DataFrame(score = div$MLModelFit.BrLen0)

hits <- findOverlaps(dm.gr, div.gr) 
ranges(div.gr)[subjectHits(hits)] = ranges(dm.gr)[queryHits(hits)]

lands.gr.df <- as.data.frame(dm.gr)
div.gr.df <- as.data.frame(div.gr)
# for some reason GRanges has problems converting start coordinates within 6 lands bins
div.gr.df <- div.gr.df[which(((div.gr.df$width - 1) %% 200000) == 0),]

tmp <- ddply(.data = div.gr.df[-c(1,5)], .variables = "start", .fun = colMeans, .progress = "text")
div.gr.df.2 <- div.gr.df[!duplicated(div.gr.df$start),]
div.gr.df.2$score <- tmp$score

lands.gr.df.2 <- lands.gr.df[which(lands.gr.df$start %in% div.gr.df.2$start),]

lands.div.dm <- cbind.data.frame(lands.gr.df.2[,-(which(names(lands.gr.df.2) == "strand"))], div.gr.df.2$score)
names(lands.div.dm)[ncol(lands.div.dm)] <- "divergence"

cor.test(lands.div.dm$theta, lands.div.dm$divergence, method = "spearman")
# S = 30874, p-value = 0.0407
# sample estimates:
#  rho 
# 0.2589766 

# NOTE cannot build a linear model due to extreme multi-collinearity (huge VIFs -- not shown)
cor.mat <- cor(lands.div.dm[,5:9], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, #cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# Moving on to partial correlations...

# NS
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$tmrca, z = cbind(lands.div.dm$theta, lands.div.dm$rho), method = "spearman")
# NS
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$rho, z = cbind(lands.div.dm$theta, lands.div.dm$tmrca), method = "spearman")


# afaik there are two hypothesis for the correlation between divergence and diversity:
# first, mutation rate variation
# second, background selection
cor.test(x = lands.div.dm$divergence, y = lands.div.dm$diversity, method = "spearman")
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$diversity, lands.div.dm$theta, method = "spearman")


# 0.25, p = 0.05
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$theta, z = cbind(lands.div.dm$rho, lands.div.dm$tmrca), method = "spearman")

# plot
theta.div.scatter <- ggplot(data = lands.div.dm, aes(x = lands.div.dm$theta, y = lands.div.dm$divergence, color = lands.div.dm$tmrca))
theta.div.scatter <- theta.div.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
#theta.div.scatter <- theta.div.scatter + geom_smooth(method = "loess", colour = 'magenta', linetype = 2, se = F)
theta.div.scatter <- theta.div.scatter + theme(plot.title = element_text(size = 20),
                                               axis.title.x = element_text(size = 20),
                                               axis.title.y = element_text(size = 20))
theta.div.scatter <- theta.div.scatter + labs(title = NULL, #"Mutation rate vs Divergence D. melanogaster : D. yakuba",
                                              color = "TMRCA",
                                              x = "Scaled Mutation Rate in D. melanogaster",
                                              y = "Divergence")
theta.div.scatter <- theta.div.scatter + annotate("text", x = 0.008, y = 0.075, label = "Spearman's rank partial correlation = 0.25, p = 0.05")

ggsave("divergence.dm.200kb.pdf", plot = theta.div.scatter, dpi = 2000, width = 7, height = 7)



########################################
#
# Divergence --- 1 Mb (TOO FEW WINDOWS HERE)
#
########################################

# to get the ranges
diversity.dm <- diversity.dm.1Mb[which(intersect.1Mb == F),] 
dm.1Mb <- as.data.frame(cbind(diversity.dm[,(1:3)], dm.lands.1Mb[,1:4]))
names(dm.1Mb) <- c("Chr", "Start", "Stop", names(dm.lands.1Mb)[1:4])
dm.1Mb$Chr <- "2L"

# converts objects
dm.gr <- makeGRangesFromDataFrame(dm.1Mb)
values(dm.gr) <- dm.1Mb[,(4:7)]
div.gr <- makeGRangesFromDataFrame(div)
values(div.gr) <- DataFrame(score = div$MLModelFit.BrLen0)

hits <- findOverlaps(dm.gr, div.gr) 
ranges(div.gr)[subjectHits(hits)] = ranges(dm.gr)[queryHits(hits)]

lands.gr.df <- as.data.frame(dm.gr)
div.gr.df <- as.data.frame(div.gr)
# for some reason GRanges has problems converting start coordinates within 6 lands bins
div.gr.df <- div.gr.df[which(((div.gr.df$width - 1) %% 1000000) == 0),]

tmp <- ddply(.data = div.gr.df[-c(1,5)], .variables = "start", .fun = colMeans, .progress = "text")
div.gr.df.2 <- div.gr.df[!duplicated(div.gr.df$start),]
div.gr.df.2$score <- tmp$score

lands.gr.df.2 <- lands.gr.df[which(lands.gr.df$start %in% div.gr.df.2$start),]

lands.div.dm <- cbind.data.frame(lands.gr.df.2[,-(which(names(lands.gr.df.2) == "strand"))], div.gr.df.2$score)
names(lands.div.dm)[ncol(lands.div.dm)] <- "divergence"

cor.test(lands.div.dm$theta, lands.div.dm$divergence, method = "spearman")
# S = 30874, p-value = 0.3
#sample estimates:
#  rho 
# 0.30

# Moving on to partial correlations...

# NS
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$tmrca, z = cbind(lands.div.dm$theta, lands.div.dm$rho), method = "spearman")
# NS
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$rho, z = cbind(lands.div.dm$theta, lands.div.dm$tmrca), method = "spearman")

# 0.57, p = 0.05
pcor.test(x = lands.div.dm$divergence, y = lands.div.dm$theta, z = cbind(lands.div.dm$rho, lands.div.dm$tmrca), method = "spearman")

theta.div.scatter <- ggplot(data = lands.div.dm, aes(x = lands.div.dm$theta, y = lands.div.dm$divergence, color = lands.div.dm$tmrca))
theta.div.scatter <- theta.div.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
#theta.div.scatter <- theta.div.scatter + geom_smooth(method = "loess", colour = 'magenta', linetype = 2, se = F)
theta.div.scatter <- theta.div.scatter + theme(plot.title = element_text(size = 20),
                                               axis.title.x = element_text(size = 20),
                                               axis.title.y = element_text(size = 20))
theta.div.scatter <- theta.div.scatter + labs(title = NULL, #"Mutation rate vs Divergence D. melanogaster : D. yakuba",
                                              color = "TMRCA",
                                              x = "Scaled Mutation Rate in D. melanogaster",
                                              y = "Divergence")
theta.div.scatter <- theta.div.scatter + annotate("text", x = 0.009, y = 0.069, label = "Spearman's rank partial correlation = 0.57, p = 0.05")

ggsave("divergence.dm.1Mb.pdf", plot = theta.div.scatter, dpi = 2000, width = 7, height = 7)

########################################
#
# TMRCA in coding vs non-coding regions
#
########################################

# TODO bin TMRCA landscapes in 1KB windows for this analysis!
  
dm.genes.coord <- read.table("~/Data/iSMC/theta_paper/real_data/droso.misc/Dmel_trans_nooverlaps-final.txt", header = F)
dm.genes.coord <- dm.genes.coord[which(dm.genes.coord$V1 == "2L"), c(1, 2, 3, 5, 6)]
names(dm.genes.coord) <- c("chrom", "chromStart", "chromEnd", "gene", "length")

dm.maps.50k <- as.data.frame(cbind(diversity.dm.50kb$chrom, diversity.dm.50kb$chromStart, diversity.dm.50kb$chromEnd,
                                   diversity.dm.50kb$avg, theta.dm.50kb$sample_mean, rho.dm.50kb$sample_mean, tmrca.dm.50kb$sample_mean))
names(dm.maps.50k) <- c("chrom", "chromStart", "chromEnd", "diversity", "theta", "rho", "tmrca")
dm.maps.50k$chrom <- "2L"

# filters based on missing data ( > 50% per window)
dm.maps.50k <- dm.maps.50k[which(intersect.50kb == F),]

# grouping per window of constant rates
dm.maps.gr <- makeGRangesFromDataFrame(dm.maps.50k)
values(dm.maps.gr) <- dm.maps.50k[,(4:7)]
genes.gr <- makeGRangesFromDataFrame(dm.genes.coord)
values(genes.gr) <- dm.genes.coord[,(5)]

hits <- findOverlaps(dm.maps.gr, genes.gr) 
ranges(genes.gr)[subjectHits(hits)] <- ranges(dm.maps.gr)[queryHits(hits)]

genes.gr.df <- as.data.frame(genes.gr)
dm.maps.gr.df <- as.data.frame(dm.maps.gr)

# gets mean of a statistic for all genes in the window
tmp <- ddply(.data = genes.gr.df[-c(1,5)], .variables = "X", .fun = colMeans, .progress = "text")
genes.gr.df.2 <- tmp %>% arrange(start)

# due to previously filtering out windows based on missing data,
# some genes do not match to a window
genes.gr.df.2 <- genes.gr.df.2[which(((genes.gr.df.2$width - 1) %% 50000) == 0),]

dm.maps.gr.df.overlap <- dm.maps.gr.df[which(dm.maps.gr.df$start %in% genes.gr.df.2$start),]

`%notin%` <- Negate(`%in%`)
dm.maps.gr.df.NOT.overlap <- dm.maps.gr.df[which(dm.maps.gr.df$start %notin% genes.gr.df.2$start),]

wilcox.test(dm.maps.gr.df.overlap$tmrca, dm.maps.gr.df.NOT.overlap$tmrca)

dm.maps.genes <- merge(dm.maps.gr.df.overlap, genes.gr.df.2, by = "start")

pcor.test(dm.maps.genes$tmrca, dm.maps.genes$X, dm.maps.genes$theta, method = "spearman")
plot(x=dm.maps.genes$tmrca, y=dm.maps.genes$X)
boxplot(dm.maps.gr.df.overlap$tmrca, dm.maps.gr.df.NOT.overlap$tmrca)



########################################
#
# Mutation x Evolutionary (Protein) Rates --- 50 kb
#
########################################

# loads 
dm.raw <- read.table("~/Data/iSMC/theta_paper/real_data/droso.misc/dpgp3_Dyak_bpp.all.csv", header = T, fill = T, stringsAsFactors = T)
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

dm.genes.coord <- read.table("~/Data/iSMC/theta_paper/real_data/droso.misc/Dmel_trans_nooverlaps-final.txt", header = F)
names(dm.genes.coord) <- c("chr", "start", "end", "x", "geneID", "length")
dm.genes.coord <- dm.genes.coord[which(dm.genes.coord$chr == "2L"), c(1, 2, 3, 5, 6)]

dm.evol <- merge(dm.genes.coord, dm.tbl.popstats.clean, by = "geneID")
dm.evol <- dm.evol[order(dm.evol$start),]


# to get the ranges
diversity.dm <- diversity.dm.50kb[which(intersect.50kb == F),]
dm.50kb <- as.data.frame(cbind(diversity.dm[,(1:3)], dm.lands.50kb)) 
names(dm.50kb) <- c("chr", "start", "end", names(dm.lands.50kb))
dm.50kb$chr <- "2L"


# grouping per gene coordinate
dm.lands.gr <- makeGRangesFromDataFrame(dm.50kb)
values(dm.lands.gr) <- dm.50kb[,(4:7)]
evolrate.gr <- makeGRangesFromDataFrame(dm.evol)
values(evolrate.gr) <- dm.evol[,(5:11)]

hits <- findOverlaps(evolrate.gr, dm.lands.gr) 

evolrate.gr.df <- as.data.frame(evolrate.gr[queryHits(hits)], row.names = NULL)
dm.lands.gr.df <- as.data.frame(dm.lands.gr[subjectHits(hits)], row.names = NULL)

dm.lands.evolrate <- cbind.data.frame(dm.lands.gr.df[,c(2,3,6:9)], evolrate.gr.df[,c(2,3,6:12)])
dm.lands.evolrate <- dm.lands.evolrate[which(dm.lands.evolrate$PiNPiS < 1),]
names(dm.lands.evolrate)[1] <- "start.window"
names(dm.lands.evolrate)[2] <- "end.window"
names(dm.lands.evolrate)[7] <- "start.gene"
names(dm.lands.evolrate)[8] <- "end.gene"

# linear model in coding regions
m.diversity.cds <- lm(diversity ~ theta + rho + tmrca, data = dm.lands.evolrate)
plot(m.diversity.cds, which = 2) 
shapiro.test(resid(m.diversity.cds)) # ***
hmctest(m.diversity.cds, nsim = 3000) # NS
dwtest(m.diversity.cds) # ***

bc.diversity.cds <- boxcox(m.diversity.cds, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity.cds$x[which.max(bc.diversity.cds$y)]
m.diversity.cds.bc <- update(m.diversity.cds, (diversity^l -1)/l~.)
plot(m.diversity.cds.bc, which = 2)
shapiro.test(resid(m.diversity.cds.bc)) # *** and worse than before => discard BC transform

dm.lands.evolrate$bin <- 1:nrow(dm.lands.evolrate)
# Because of auto-correlation we compute p-values for the variables using a GLS
g.diversity <- gls(diversity ~ theta + rho + tmrca, data = dm.lands.evolrate, corr = corAR1(0, ~bin))

summary(m.diversity.cds)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.0100127  0.0002685 -37.290   <2e-16 ***
# theta        0.9371188  0.0118775  78.899   <2e-16 ***
# rho          0.0061516  0.0028670   2.146   0.0326 *  
# tmrca        0.0111534  0.0003803  29.330   <2e-16 ***

summary(g.diversity)
#  Phi 
# 0.7088267 

# Coefficients:
#   Value   Std.Error   t-value p-value
# (Intercept) -0.0111820 0.000267317 -41.83045  0.0000
# theta        0.9457321 0.016595598  56.98693  0.0000
# rho          0.0057487 0.003302340   1.74080  0.0827
# tmrca        0.0123129 0.000370546  33.22915  0.0000

# type 2 anova
anova.diversity.cds <- Anova(m.diversity.cds)
apiss <- anova.diversity.cds$"Sum Sq"
anova.diversity.cds$VarExp <- apiss / sum(apiss)

anova.diversity.cds
# Anova Table (Type II tests)

# Response: diversity
# Sum Sq  Df   F value   Pr(>F)  VarExp
# theta     0.00077936   1 6224.9834 0.000000 0.83896
# rho       0.00000058   1    4.6039 0.032629 0.00062
# tmrca     0.00010770   1  860.2648 0.000000 0.11594
# Residuals 0.00004132 330                    0.04448

#
cor.test(dm.lands.evolrate$PiNPiS, dm.lands.evolrate$rho, method = "spearman")
pcor.test(dm.lands.evolrate$PiNPiS, dm.lands.evolrate$theta, dm.lands.evolrate$rho, method = "spearman") 
pcor.test(dm.lands.evolrate$PiNPiS, dm.lands.evolrate$theta, dm.lands.evolrate$tmrca, method = "spearman") 

# decomposing PiN/PiS
pcor.test(dm.lands.evolrate$PiN, dm.lands.evolrate$theta, dm.lands.evolrate$tmrca, method = "spearman") 
pcor.test(dm.lands.evolrate$PiS, dm.lands.evolrate$theta, dm.lands.evolrate$tmrca, method = "spearman") 

cor.test(dm.lands.evolrate$PiS, dm.lands.evolrate$rho, method = "spearman") 
cor.test(dm.lands.evolrate$dS, dm.lands.evolrate$rho, method = "spearman") 

cor.test(dm.lands.evolrate$dS, dm.lands.evolrate$PiS, method = "spearman")
pcor.test(dm.lands.evolrate$dS, dm.lands.evolrate$theta, dm.lands.evolrate$tmrca, method = "spearman")
pcor.test(dm.lands.evolrate$dS, dm.lands.evolrate$PiS, dm.lands.evolrate$theta, method = "spearman")

mean(dm.lands.evolrate$dS)
mean(div$MLModelFit.BrLen0)

sd(dm.lands.evolrate$theta) / mean(dm.lands.evolrate$theta)
sd(dm.50kb$theta) / mean(dm.50kb$theta)

#####################
#
# B-value statistic across the Drosophila Genome
#
#####################

# TODO


