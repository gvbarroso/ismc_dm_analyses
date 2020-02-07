# Created: 24/06/2019
# Last modified: 28/06/2019
# Author: Gustavo Barroso

library(ppcor)
library(MASS)
library(reshape2)
library(ggplot2)
library(ppls)
library(scales)
library(cowplot)
library(MASS)
library(lmtest)
library(relaimpo)
library(corrplot)
library(nlme)
library(car)
library(GenomicRanges)


croll <- read.table("Data/iSMC/rho_paper/real_data/3_species_maps/CrollEtAl/cross_maps_20kb.txt")

# rho-modulated 40x5x1
mono.rho.20kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/40x5x1/rho.zt.rho.joint.20kb.120000-5999999.bedgraph", header = T)
cor.test(x = croll$cross_avg, y = mono.rho.20kb$sample_mean, method = "pearson")$estimate ^ 2

# R_2 table for plotting at the end
r2.tab <- as.data.frame(matrix(ncol = 5, nrow = 3))
colnames(r2.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Bin_size(kb)")

########################################
#
# 30x5x5 --- 20kb
#
########################################

# 20kb 

# recombination landscapes
rho.zt.20kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.rho.20kb.bedgraph", header = T)
# R2 with genetic map increases relative to 40x5x1 model with same sample size
cor.test(x = croll$cross_avg, y = rho.zt.20kb$sample_mean, method = "pearson")$estimate ^ 2

# diversity
diversity.zt.20kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.diversity.20kb.bedgraph", header = T)
diversity.zt.20kb$avg <- apply(diversity.zt.20kb[4:ncol(diversity.zt.20kb)], 1, mean)

# mutation landscapes
theta.zt.20kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.theta.20kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.zt.20kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.TMRCA.20kb.bedgraph", header = T)
tmrca.zt.20kb$sample_mean <- apply(tmrca.zt.20kb[4:ncol(tmrca.zt.20kb)], 1, mean)

# missing data
missing.prop.20kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.missing.prop.20kb.bedgraph", header = T)
intersect.20kb <- apply(missing.prop.20kb[4:ncol(missing.prop.20kb)], 1, function(x) any(x > 0.25)) 

zt.lands.20kb <- as.data.frame(cbind(diversity.zt.20kb$avg, theta.zt.20kb$sample_mean, rho.zt.20kb$sample_mean, tmrca.zt.20kb$sample_mean))
names(zt.lands.20kb) <- c("diversity", "theta", "rho", "tmrca")
zt.lands.20kb$bin <- 1:nrow(zt.lands.20kb)
# filters based on missing data ( > 50% per window)
zt.lands.20kb <- zt.lands.20kb[which(intersect.20kb == F),]

# multi-collinearity, unlike neutral simulations
cor.mat <- cor(zt.lands.20kb[,1:4], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", tl.cex = 1.5, is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# we try to disentangle using partial correlations
pcor.mat <- cor.mat

pcor.test(x = zt.lands.20kb$rho, y = zt.lands.20kb$tmrca, z = zt.lands.20kb$theta, method = "spearman")
pcor.mat[4, 3] <- 0.48
pcor.test(x = zt.lands.20kb$rho, y = zt.lands.20kb$diversity, z = zt.lands.20kb$tmrca, method = "spearman")
pcor.mat[3, 1] <- 0
pcor.test(x = zt.lands.20kb$rho, y = zt.lands.20kb$theta, z = zt.lands.20kb$tmrca, method = "spearman")
pcor.mat[3, 2] <- 0

pcor.test(x = zt.lands.20kb$theta, y = zt.lands.20kb$diversity, z = zt.lands.20kb$tmrca, method = "spearman")
pcor.test(x = zt.lands.20kb$theta, y = zt.lands.20kb$diversity, z = zt.lands.20kb$rho, method = "spearman")
pcor.test(x = zt.lands.20kb$theta, y = zt.lands.20kb$diversity, z = cbind(zt.lands.20kb$rho, zt.lands.20kb$tmrca), method = "spearman")
pcor.mat[2, 1] <- 0.69
pcor.test(x = zt.lands.20kb$theta, y = zt.lands.20kb$tmrca, z = zt.lands.20kb$rho, method = "spearman")
pcor.mat[4, 2] <- 0.45

pcor.test(x = zt.lands.20kb$tmrca, y = zt.lands.20kb$diversity, z = zt.lands.20kb$rho, method = "spearman")
pcor.test(x = zt.lands.20kb$tmrca, y = zt.lands.20kb$diversity, z = zt.lands.20kb$theta, method = "spearman")
pcor.test(x = zt.lands.20kb$tmrca, y = zt.lands.20kb$diversity, z = cbind(zt.lands.20kb$rho, zt.lands.20kb$theta), method = "spearman")
pcor.mat[4, 1] <- 0.96

pdf("Data/iSMC/theta_paper/submission/Figures/zt.pcor.20kb.pdf")
corrplot(pcor.mat, method = "color", addCoef.col = "black", tl.cex = 1.5, is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()

# CV (added to plots with annotate)
sd(zt.lands.20kb$diversity) / mean(zt.lands.20kb$diversity) # 0.74
sd(zt.lands.20kb$theta) / mean(zt.lands.20kb$theta) # 0.09
sd(zt.lands.20kb$rho) / mean(zt.lands.20kb$rho) # 0.31
sd(zt.lands.20kb$tmrca) / mean(zt.lands.20kb$tmrca) # 0.53

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(zt.lands.20kb[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 20, y = value)) 
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
diversity.map <- diversity.map + annotate("text", x = 3000, y = 0.055, label = "CV = 74%")

molten.rho <- melt(zt.lands.20kb[c(3,5)], id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 20, y = value)) 
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 3000, y = 0.05, label = "CV = 31%")

molten.theta <- melt(zt.lands.20kb[c(2,5)], id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 20, y = value)) 
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 3000, y = 0.019, label = "CV = 9%")

molten.tmrca <- melt(zt.lands.20kb[c(4,5)], id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 20, y = value)) 
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
tmrca.map <- tmrca.map + annotate("text", x = 3000, y = 2.8, label = "CV = 53%")

p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
p


cowplot::ggsave("Zt.maps.20kb.pdf", plot = p, device = "pdf", dpi = 500, width = 10, height = 8,
                path = "Data/iSMC/theta_paper/submission/Figures/")





# Linear model with all variables
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = zt.lands.20kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) #***
hist(resid(m.diversity.bc), nclass = 30) # not too bad?
hmctest(m.diversity.bc, nsim = 3000) # *
dwtest(m.diversity.bc) # NS
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# (Intercept) -1.594071   0.007045 -226.259   <2e-16 ***
# theta        7.194686   0.571384   12.592   <2e-16 ***
# rho          0.137791   0.087334    1.578    0.116    
# tmrca        0.057717   0.001277   45.196   <2e-16 ***
# Multiple R-squared:  0.953,	Adjusted R-squared:  0.9524 

# type 2 anova
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# theta     0.013263   1  158.5505 0.00000 0.06493
# rho       0.000208   1    2.4893 0.11595 0.00102
# tmrca     0.170872   1 2042.6638 0.00000 0.83657
# Residuals 0.019909 238                   0.09747

r2.tab[1,] <- c(6.5 + 83.7, 6.5, 0, 83.7, 20)









# Linear model without tmrca -> rho becomes significant
m.diversity <- lm(diversity ~ rho + theta, data = zt.lands.20kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.0, 1.0, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # **
hist(resid(m.diversity.bc), nclass = 30) # not too bad?
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # ***
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# (Intercept) -2.91243    0.06025 -48.336  < 2e-16 ***
# rho          4.73920    0.79872   5.933 1.03e-08 ***
# theta       53.45874    4.72471  11.315  < 2e-16 ***
# Multiple R-squared:  0.5239,	Adjusted R-squared:  0.5199 

anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# rho         1 1.02188 1.02188  134.97 4.9530e-25 0.26887
# theta       1 0.96927 0.96927  128.02 4.7366e-24 0.25503
# Residuals 239 1.80950 0.00757                    0.47610


########################################
#
# 30x5x5 --- 50kb
#
########################################

# 50kb 

# recombination landscapes
rho.zt.50kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.rho.50kb.bedgraph", header = T)

# diversity
diversity.zt.50kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.diversity.50kb.bedgraph", header = T)
diversity.zt.50kb$avg <- apply(diversity.zt.50kb[4:ncol(diversity.zt.50kb)], 1, mean)

# mutation landscapes
theta.zt.50kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.theta.50kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.zt.50kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.TMRCA.50kb.bedgraph", header = T)
tmrca.zt.50kb$sample_mean <- apply(tmrca.zt.50kb[4:ncol(tmrca.zt.50kb)], 1, mean)

# missing data
missing.prop.50kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.missing.prop.50kb.bedgraph", header = T)
intersect.50kb <- apply(missing.prop.50kb[4:ncol(missing.prop.50kb)], 1, function(x) any(x > 0.25)) 

zt.lands.50kb <- as.data.frame(cbind(diversity.zt.50kb$avg, theta.zt.50kb$sample_mean, rho.zt.50kb$sample_mean, tmrca.zt.50kb$sample_mean))
names(zt.lands.50kb) <- c("diversity", "theta", "rho", "tmrca")
zt.lands.50kb$bin <- 1:nrow(zt.lands.50kb)
# filters based on missing data ( > 50% per window)
zt.lands.50kb <- zt.lands.50kb[which(intersect.50kb == F),]

# multi-collinearity, unlike neutral simulations
cor.mat <- cor(zt.lands.50kb[,1:4], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = .5, diag = F, type = "lower")

# we try to disentangle using partial correlations
pcor.mat <- cor.mat

pcor.test(x = zt.lands.50kb$rho, y = zt.lands.50kb$tmrca, z = zt.lands.50kb$theta, method = "spearman")
pcor.mat[4, 3] <- 0.43
pcor.test(x = zt.lands.50kb$rho, y = zt.lands.50kb$diversity, z = zt.lands.50kb$tmrca, method = "spearman")
pcor.mat[3, 1] <- 0
pcor.test(x = zt.lands.50kb$rho, y = zt.lands.50kb$theta, z = zt.lands.50kb$tmrca, method = "spearman")
pcor.mat[3, 2] <- 0

pcor.test(x = zt.lands.50kb$theta, y = zt.lands.50kb$diversity, z = zt.lands.50kb$tmrca, method = "spearman")
pcor.test(x = zt.lands.50kb$theta, y = zt.lands.50kb$diversity, z = zt.lands.50kb$rho, method = "spearman")
pcor.test(x = zt.lands.50kb$theta, y = zt.lands.50kb$diversity, z = cbind(zt.lands.50kb$rho, zt.lands.50kb$tmrca), method = "spearman")
pcor.mat[2, 1] <- 0.66
pcor.test(x = zt.lands.50kb$theta, y = zt.lands.50kb$tmrca, z = zt.lands.50kb$rho, method = "spearman")
pcor.mat[4, 2] <- 0.56

pcor.test(x = zt.lands.50kb$tmrca, y = zt.lands.50kb$diversity, z = zt.lands.50kb$rho, method = "spearman")
pcor.test(x = zt.lands.50kb$tmrca, y = zt.lands.50kb$diversity, z = zt.lands.50kb$theta, method = "spearman")
pcor.test(x = zt.lands.50kb$tmrca, y = zt.lands.50kb$diversity, z = cbind(zt.lands.50kb$theta, zt.lands.50kb$rho), method = "spearman")
pcor.mat[4, 1] <- 0.94

pdf("Data/iSMC/theta_paper/submission/Figures/zt.pcor.50kb.pdf")
corrplot(pcor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()

# CV (added to plots with annotate)
sd(zt.lands.50kb$diversity) / mean(zt.lands.50kb$diversity) # 0.53
sd(zt.lands.50kb$theta) / mean(zt.lands.50kb$theta) # 0.06
sd(zt.lands.50kb$rho) / mean(zt.lands.50kb$rho) # 0.19
sd(zt.lands.50kb$tmrca) / mean(zt.lands.50kb$tmrca) # 0.38

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(zt.lands.50kb[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 50, y = value)) 
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
diversity.map <- diversity.map + annotate("text", x = 3000, y = 0.035, label = "CV = 53%")

molten.rho <- melt(zt.lands.50kb[c(3,5)], id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 50, y = value)) 
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 50), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 3000, y = 0.035, label = "CV = 19%")

molten.theta <- melt(zt.lands.50kb[c(2,5)], id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 50, y = value)) 
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 3000, y = 0.016, label = "CV = 6%")

molten.tmrca <- melt(zt.lands.50kb[c(4,5)], id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 50, y = value)) 
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
tmrca.map <- tmrca.map + annotate("text", x = 3000, y = 2.25, label = "CV = 38%")

p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
p


cowplot::ggsave("Zt.maps.50kb.pdf", plot = p, device = "pdf", dpi = 500, width = 10, height = 8,
                path = "Data/iSMC/theta_paper/submission/Figures/")





# Linear model with all variables
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = zt.lands.50kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.2, 1.2, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) #***
hist(resid(m.diversity.bc), nclass = 30) # not too bad?
hmctest(m.diversity.bc, nsim = 3000) # *
dwtest(m.diversity.bc) # NS
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# (Intercept) -1.248433   0.006102 -204.591  < 2e-16 ***
# theta        4.158371   0.492833    8.438 3.13e-13 ***
# rho          0.060573   0.078012    0.776    0.439    
# tmrca        0.027706   0.001023   27.092  < 2e-16 ***
# Multiple R-squared:  0.9548,	Adjusted R-squared:  0.9534 

# type 2 anova
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# theta     0.0006996  1  71.1946 0.00000 0.07886
# rho       0.0000059  1   0.6029 0.43937 0.00067
# tmrca     0.0072123  1 733.9644 0.00000 0.81302
# Residuals 0.0009532 97                  0.10745

r2.tab[2,] <- c(81.3 + 7.9, 7.9, 0, 81.3, 50)








# Linear model without tmrca
m.diversity <- lm(diversity ~ rho + theta, data = zt.lands.50kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.0, 1.0, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # NS
hist(resid(m.diversity.bc), nclass = 30) 
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # *
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# rho        1 0.17938 0.179385  74.046 1.2769e-13 0.29912
# theta      1 0.18291 0.182908  75.501 8.4090e-14 0.30500
# Residuals 98 0.23741 0.002423                    0.39588



########################################
#
# 30x5x5 --- 200kb
#
########################################

# recombination landscapes
rho.zt.200kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.rho.200kb.bedgraph", header = T)

# diversity
diversity.zt.200kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.diversity.200kb.bedgraph", header = T)
diversity.zt.200kb$avg <- apply(diversity.zt.200kb[4:ncol(diversity.zt.200kb)], 1, mean)

# mutation landscapes
theta.zt.200kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.theta.200kb.bedgraph", header = T)

# TMRCA landscapes
tmrca.zt.200kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.TMRCA.200kb.bedgraph", header = T)
tmrca.zt.200kb$sample_mean <- apply(tmrca.zt.200kb[4:ncol(tmrca.zt.200kb)], 1, mean)

# missing data
missing.prop.200kb <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/30x5x5/zt.30x5x5.missing.prop.200kb.bedgraph", header = T)
intersect.200kb <- apply(missing.prop.200kb[4:ncol(missing.prop.200kb)], 1, function(x) any(x > 0.25)) 

zt.lands.200kb <- as.data.frame(cbind(diversity.zt.200kb$avg, theta.zt.200kb$sample_mean, rho.zt.200kb$sample_mean, tmrca.zt.200kb$sample_mean))
names(zt.lands.200kb) <- c("diversity", "theta", "rho", "tmrca")
zt.lands.200kb$bin <- 1:nrow(zt.lands.200kb)
# filters based on missing data ( > 50% per window)
zt.lands.200kb <- zt.lands.200kb[which(intersect.200kb == F),]

# multi-collinearity, unlike neutral simulations
cor.mat <- cor(zt.lands.200kb[,1:4], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# we try to disentangle using partial correlations
pcor.mat <- cor.mat
  
pcor.test(x = zt.lands.200kb$rho, y = zt.lands.200kb$tmrca, z = zt.lands.200kb$theta, method = "spearman")
pcor.mat[4, 3] <- 0.57
pcor.test(x = zt.lands.200kb$rho, y = zt.lands.200kb$diversity, z = zt.lands.200kb$tmrca, method = "spearman")
pcor.mat[3, 1] <- 0
pcor.test(x = zt.lands.200kb$rho, y = zt.lands.200kb$theta, z = zt.lands.200kb$tmrca, method = "spearman")
pcor.mat[3, 2] <- 0

pcor.test(x = zt.lands.200kb$theta, y = zt.lands.200kb$diversity, z = zt.lands.200kb$tmrca, method = "spearman")
pcor.test(x = zt.lands.200kb$theta, y = zt.lands.200kb$diversity, z = zt.lands.200kb$rho, method = "spearman")
pcor.test(x = zt.lands.200kb$theta, y = zt.lands.200kb$diversity, z = cbind(zt.lands.200kb$tmrca, zt.lands.200kb$rho), method = "spearman")
pcor.mat[2, 1] <- 0.73
pcor.test(x = zt.lands.200kb$theta, y = zt.lands.200kb$tmrca, z = zt.lands.200kb$rho, method = "spearman")
pcor.mat[4, 2] <- 0.56

pcor.test(x = zt.lands.200kb$tmrca, y = zt.lands.200kb$diversity, z = zt.lands.200kb$rho, method = "spearman")
pcor.test(x = zt.lands.200kb$tmrca, y = zt.lands.200kb$diversity, z = zt.lands.200kb$theta, method = "spearman")
pcor.test(x = zt.lands.200kb$tmrca, y = zt.lands.200kb$diversity, z = cbind(zt.lands.200kb$theta, zt.lands.200kb$rho), method = "spearman")
pcor.mat[4, 1] <- 0.88

pdf("Data/iSMC/theta_paper/submission/Figures/zt.pcor.200kb.pdf")
corrplot(pcor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")
dev.off()

# CV (added to plots with annotate)
sd(zt.lands.200kb$diversity) / mean(zt.lands.200kb$diversity) # 0.40
sd(zt.lands.200kb$theta) / mean(zt.lands.200kb$theta) # 0.05
sd(zt.lands.200kb$rho) / mean(zt.lands.200kb$rho) # 0.13
sd(zt.lands.200kb$tmrca) / mean(zt.lands.200kb$tmrca) # 0.28

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)

molten.diversity <- melt(zt.lands.200kb[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 200, y = value)) 
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
diversity.map <- diversity.map + annotate("text", x = 3000, y = 0.026, label = "CV = 40%")

molten.rho <- melt(zt.lands.200kb[c(3,5)], id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 200, y = value)) 
rho.map <- rho.map + geom_line(data = molten.rho, colour = "#7CAE00")
rho.map <- rho.map + geom_smooth(method = "loess", se = F, colour = "#7CAE00")
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 3000, y = 0.032, label = "CV = 13%")

molten.theta <- melt(zt.lands.200kb[c(2,5)], id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 200, y = value)) 
theta.map <- theta.map + geom_line(data = molten.theta, colour = "#00BFC4")
theta.map <- theta.map + geom_smooth(method = "loess", se = F, colour = "#00BFC4")
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta)) 
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 3000, y = 0.0145, label = "CV = 5%")

molten.tmrca <- melt(zt.lands.200kb[c(4,5)], id.vars = "bin")
tmrca.map <- ggplot(data = molten.tmrca, aes(x = bin * 200, y = value)) 
tmrca.map <- tmrca.map + geom_line(data = molten.tmrca, colour = "#88419d")
tmrca.map <- tmrca.map + geom_smooth(method = "loess", se = F, colour = "#88419d")
tmrca.map <- tmrca.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
tmrca.map <- tmrca.map + labs(title = NULL, x = "Position (kb)", y = expression(tau))
tmrca.map <- tmrca.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
tmrca.map <- tmrca.map + annotate("text", x = 3000, y = 1.7, label = "CV = 28%")

p <- plot_grid(diversity.map, theta.map, rho.map, tmrca.map, nrow = 4, ncol = 1, labels = NULL, label_size = 18, scale = 0.9)
p

cowplot::ggsave("Zt.maps.200kb.pdf", plot = p, device = "pdf", dpi = 500, width = 10, height = 8,
                path = "Data/iSMC/theta_paper/submission/Figures/")




# Linear model with all variables
m.diversity <- lm(diversity ~ theta + rho + tmrca, data = zt.lands.200kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.7, 1.7, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # NS
hist(resid(m.diversity.bc), nclass = 30) 
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # NS
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

summary(m.diversity.bc)
# (Intercept) -0.8562889  0.0027544 -310.885  < 2e-16 ***
# theta        1.0370580  0.2440874    4.249 0.000303 ***
# rho          0.0208640  0.0439954    0.474 0.639806    
# tmrca        0.0063777  0.0005374   11.869 2.75e-11 ***
# Multiple R-squared:  0.9659,	Adjusted R-squared:  0.9615

vif(m.diversity.bc)
# theta      rho    tmrca 
# 2.781900 1.890101 2.626465 

# type 2 anova
anova.diversity <- Anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# theta     4.785e-06  1  18.0516 0.00030 0.09911
# rho       6.000e-08  1   0.2249 0.63981 0.00123
# tmrca     3.734e-05  1 140.8615 0.00000 0.77338
# Residuals 6.097e-06 23                  0.12628

r2.tab[3,] <- c(77.3 + 9.9, 9.9, 0, 77.3, 200)







# Linear model without tmrca
m.diversity <- lm(diversity ~ rho + theta, data = zt.lands.200kb)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0.0, 1.0, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) # NS
hist(resid(m.diversity.bc), nclass = 30) # not too bad?
hmctest(m.diversity.bc, nsim = 3000) # NS
dwtest(m.diversity.bc) # *
Box.test(resid(m.diversity.bc)[order(predict(m.diversity.bc))], type = "Ljung-Box") # NS

anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# rho        1 0.022765 0.0227652  53.977 2.4389e-09 0.36516
# theta      1 0.019755 0.0197553  46.841 1.4095e-08 0.31688
# Residuals 47 0.019822 0.0004218                    0.31796



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
r2.plot <- r2.plot + scale_x_continuous(breaks = c(20, 50, 200)) 
r2.plot <- r2.plot + scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 100))
r2.plot <- r2.plot + labs(title = NULL, x = "Bin Size (kb)", y = "Var. Explained")
r2.plot <- r2.plot + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

ggsave("Data/iSMC/theta_paper/submission/Figures/lm.r2.zt.pdf", plot = r2.plot, dpi = 500, width = 7, height = 7)


########################################
#
# Divergence --- 20kb
#
########################################

# divergence data from Z. tritici and Z. ardabiliae
div <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/Zt09-ST04IR111_windows_realigned_cleaned_divergence.statistics.csv", header = T)
div <- div[which(div$Chr == "chr_1"), c(1:3, 6)] # keeping only chr1 and the relevant columns

# to get the ranges
diversity.zt <- diversity.zt.20kb[which(intersect.20kb == F),] 
zt.20kb <- as.data.frame(cbind(diversity.zt[,(1:3)], zt.lands.20kb[,1:4]))
names(zt.20kb) <- c("Chr", "Start", "Stop", names(zt.lands.20kb)[1:4])
zt.20kb$Chr <- "chr_1"

# converts objects
zt.gr <- makeGRangesFromDataFrame(zt.20kb)
values(zt.gr) <- zt.20kb[,(4:7)]
div.gr <- makeGRangesFromDataFrame(div)
values(div.gr) <- DataFrame(score = div$MLModelFit.BrLen0)
  
hits <- findOverlaps(zt.gr, div.gr) 
ranges(div.gr)[subjectHits(hits)] = ranges(zt.gr)[queryHits(hits)]

div.gr.df <- as.data.frame(div.gr)
lands.gr.df <- as.data.frame(zt.gr)

# for some reason GRanges has problems converting start coordinates within 6 lands bins
div.gr.df <- div.gr.df[which(((div.gr.df$width - 1) %% 20000) == 0),]

tmp <- ddply(.data = div.gr.df[-c(1,5)], .variables = "start", .fun = colMeans, .progress = "text")
div.gr.df.2 <- div.gr.df[!duplicated(div.gr.df$start),]
div.gr.df.2$score <- tmp$score

lands.gr.df.2 <- lands.gr.df[which(lands.gr.df$start %in% div.gr.df.2$start),]

lands.div.zt <- cbind.data.frame(lands.gr.df.2[,-(which(names(lands.gr.df.2) == "strand"))], div.gr.df.2$score)
names(lands.div.zt)[ncol(lands.div.zt)] <- "divergence"

cor.test(lands.div.zt$divergence, lands.div.zt$theta, method = "spearman")
# S = 179168, p-value = 0.003233
#  rho 
# 0.2743444 

# NOTE cannot build a linear model due to extreme multi-collinearity (huge VIFs -- not shown)
cor.mat <- cor(lands.div.zt[,5:9], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# Moving on to partial correlations...

# NS
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$tmrca, z = cbind(lands.div.zt$theta, lands.div.zt$rho), method = "spearman")
# 0.23, p = 0.014
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$rho, z = cbind(lands.div.zt$theta, lands.div.zt$tmrca), method = "spearman")

# afaik there are two hypothesis for the correlation between divergence and diversity:
# first, mutation rate variation
# second, background selection
cor.test(x = lands.div.zt$divergence, y = lands.div.zt$diversity, method = "spearman")
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$diversity, lands.div.zt$theta, method = "spearman")

# 0.26, p = 0.006
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$theta, z = cbind(lands.div.zt$rho, lands.div.zt$tmrca), method = "spearman")

theta.div.scatter <- ggplot(data = lands.div.zt, aes(x = lands.div.zt$theta, y = lands.div.zt$divergence, color = lands.div.zt$tmrca))
theta.div.scatter <- theta.div.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
#theta.div.scatter <- theta.div.scatter + geom_smooth(method = "loess", colour = 'magenta', linetype = 2, se = F)
theta.div.scatter <- theta.div.scatter + theme(plot.title = element_text(size = 20),
                                               axis.title.x = element_text(size = 16),
                                               axis.title.y = element_text(size = 16))
theta.div.scatter <- theta.div.scatter + labs(title = NULL,#"Mutation rate vs Divergence Z. tritici : Z. ardabiliae",
                                              color = "TMRCA",
                                              x = "Scaled Mutation Rate in Z. tritici",
                                              y = "Divergence")
theta.div.scatter <- theta.div.scatter + annotate("text", x = 0.01385, y = 0.08, label = "Spearman's rank partial correlation = 0.26, p = 0.006")

ggsave("Data/iSMC/theta_paper/submission/Figures/divergence.zt.20kb.pdf", plot = theta.div.scatter, dpi = 500, width = 7, height = 7)


########################################
#
# Divergence --- 50kb
#
########################################

# divergence data from Z. tritici and Z. ardabiliae
div <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/Zt09-ST04IR111_windows_realigned_cleaned_divergence.statistics.csv", header = T)
div <- div[which(div$Chr == "chr_1"), c(1:3, 6)] # keeping only chr1 and the relevant columns

# to get the ranges
diversity.zt <- diversity.zt.50kb[which(intersect.50kb == F),] 
zt.50kb <- as.data.frame(cbind(diversity.zt[,(1:3)], zt.lands.50kb[,1:4]))
names(zt.50kb) <- c("Chr", "Start", "Stop", names(zt.lands.50kb)[1:4])
zt.50kb$Chr <- "chr_1"

# converts objects
zt.gr <- makeGRangesFromDataFrame(zt.50kb)
values(zt.gr) <- zt.50kb[,(4:7)]
div.gr <- makeGRangesFromDataFrame(div)
values(div.gr) <- DataFrame(score = div$MLModelFit.BrLen0)

hits <- findOverlaps(zt.gr, div.gr) 
ranges(div.gr)[subjectHits(hits)] = ranges(zt.gr)[queryHits(hits)]

div.gr.df <- as.data.frame(div.gr)
lands.gr.df <- as.data.frame(zt.gr)

# for some reason GRanges has problems converting start coordinates within 6 lands bins
div.gr.df <- div.gr.df[which(((div.gr.df$width - 1) %% 50000) == 0),]

tmp <- ddply(.data = div.gr.df[-c(1,5)], .variables = "start", .fun = colMeans, .progress = "text")
div.gr.df.2 <- div.gr.df[!duplicated(div.gr.df$start),]
div.gr.df.2$score <- tmp$score

lands.gr.df.2 <- lands.gr.df[which(lands.gr.df$start %in% div.gr.df.2$start),]

lands.div.zt <- cbind.data.frame(lands.gr.df.2[,-(which(names(lands.gr.df.2) == "strand"))], div.gr.df.2$score)
names(lands.div.zt)[ncol(lands.div.zt)] <- "divergence"

cor.test(lands.div.zt$divergence, lands.div.zt$theta, method = "spearman")
# S = 44950, p-value = 0.001581
# sample estimates:
#   rho 
# 0.3605974 

# NOTE cannot build a linear model due to extreme multi-collinearity (huge VIFs -- not shown)
cor.mat <- cor(lands.div.zt[0,5:9], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# Moving on to partial correlations...

# NS
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$tmrca, z = cbind(lands.div.zt$theta, lands.div.zt$rho), method = "spearman")
# NS
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$rho, z = cbind(lands.div.zt$theta, lands.div.zt$tmrca), method = "spearman")

# afaik there are two hypothesis for the correlation between divergence and diversity:
# first, mutation rate variation
# second, background selection
cor.test(x = lands.div.zt$divergence, y = lands.div.zt$diversity, method = "spearman")
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$diversity, lands.div.zt$theta, method = "spearman")

# 0.26, p = 0.029
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$theta, z = cbind(lands.div.zt$rho, lands.div.zt$tmrca), method = "spearman")

theta.div.scatter <- ggplot(data = lands.div.zt, aes(x = lands.div.zt$theta, y = lands.div.zt$divergence, color = lands.div.zt$tmrca))
theta.div.scatter <- theta.div.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
#theta.div.scatter <- theta.div.scatter + geom_smooth(method = "loess", colour = 'magenta', linetype = 2, se = F)
theta.div.scatter <- theta.div.scatter + theme(plot.title = element_text(size = 20),
                                               axis.title.x = element_text(size = 16),
                                               axis.title.y = element_text(size = 16))
theta.div.scatter <- theta.div.scatter + labs(title = NULL,#"Mutation rate vs Divergence Z. tritici : Z. ardabiliae",
                                              color = "TMRCA",
                                              x = "Scaled Mutation Rate in Z. tritici",
                                              y = "Divergence")
theta.div.scatter <- theta.div.scatter + annotate("text", x = 0.0145, y = 0.1, label = "Spearman's rank partial correlation = 0.26, p = 0.029")

ggsave("Data/iSMC/theta_paper/submission/Figures/divergence.zt.50kb.pdf", plot = theta.div.scatter, dpi = 500, width = 7, height = 7)


########################################
#
# Divergence --- 200kb
#
########################################

# divergence data from Z. tritici and Z. ardabiliae
div <- read.table("Data/iSMC/theta_paper/real_data/zt.chr1.theta.rho/Zt09-ST04IR111_windows_realigned_cleaned_divergence.statistics.csv", header = T)
div <- div[which(div$Chr == "chr_1"), c(1:3, 6)] # keeping only chr1 and the relevant columns

# to get the ranges
diversity.zt <- diversity.zt.200kb[which(intersect.200kb == F),]
zt.200kb <- as.data.frame(cbind(diversity.zt[,(1:3)], zt.lands.200kb[,1:4])) 
names(zt.200kb) <- c("Chr", "Start", "Stop", names(zt.lands.200kb[1:4]))
zt.200kb$Chr <- "chr_1"

# converts objects
zt.gr <- makeGRangesFromDataFrame(zt.200kb)
values(zt.gr) <- zt.200kb[,(4:7)]
div.gr <- makeGRangesFromDataFrame(div)
values(div.gr) <- DataFrame(score = div$MLModelFit.BrLen0)

hits <- findOverlaps(zt.gr, div.gr) 
ranges(div.gr)[subjectHits(hits)] = ranges(zt.gr)[queryHits(hits)]

div.gr.df <- as.data.frame(div.gr)
lands.gr.df <- as.data.frame(zt.gr)

# for some reason GRanges has problems converting start coordinates within 6 theta bins
div.gr.df <- div.gr.df[which(((div.gr.df$width - 1) %% 200000) == 0),]

tmp <- ddply(.data = div.gr.df[-c(1,5)], .variables = "start", .fun = colMeans, .progress = "text")
div.gr.df.2 <- div.gr.df[!duplicated(div.gr.df$start),]
div.gr.df.2$score <- tmp$score

lands.gr.df.2 <- lands.gr.df[which(lands.gr.df$start %in% div.gr.df.2$start),]

lands.div.zt <- cbind.data.frame(lands.gr.df.2[,-(which(names(lands.gr.df.2) == "strand"))], div.gr.df.2$score)
names(lands.div.zt)[ncol(lands.div.zt)] <- "divergence"

cor.test(lands.div.zt$divergence, lands.div.zt$theta, method = "spearman")
# S = 710, p-value = 1.317e-05
# sample estimates:
#  rho 
# 0.757265 

# NOTE cannot build a linear model due to extreme multi-collinearity (huge VIFs -- not shown)
cor.mat <- cor(lands.div.zt[,5:9], method = "spearman")
corrplot(cor.mat, method = "color", addCoef.col = "black", is.corr = T, cl.lim = c(0, 1),
         tl.col = "black", number.cex = 1.5, diag = F, type = "lower")

# Moving on to partial correlations...

# NS
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$tmrca, z = cbind(lands.div.zt$theta, lands.div.zt$rho), method = "spearman")
# NS
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$rho, z = cbind(lands.div.zt$theta, lands.div.zt$tmrca), method = "spearman")

# afaik there are two hypothesis for the correlation between divergence and diversity:
# first, mutation rate variation
# second, background selection
cor.test(x = lands.div.zt$divergence, y = lands.div.zt$diversity, method = "spearman")
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$diversity, lands.div.zt$theta, method = "spearman")

# 0.64, p = 0.0008
pcor.test(x = lands.div.zt$divergence, y = lands.div.zt$theta, z = cbind(lands.div.zt$rho, lands.div.zt$tmrca), method = "spearman")


theta.div.scatter <- ggplot(data = lands.div.zt, aes(x = lands.div.zt$theta, y = lands.div.zt$divergence, color = lands.div.zt$tmrca))
theta.div.scatter <- theta.div.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
#theta.div.scatter <- theta.div.scatter + geom_smooth(method = "loess", colour = 'magenta', linetype = 2, se = F)
theta.div.scatter <- theta.div.scatter + theme(plot.title = element_text(size = 20),
                                               axis.title.x = element_text(size = 16),
                                               axis.title.y = element_text(size = 16))
theta.div.scatter <- theta.div.scatter + labs(title = NULL,#"Mutation rate vs Divergence Z. tritici : Z. ardabiliae",
                                              color = "TMRCA",
                                              x = "Scaled Mutation Rate in Z. tritici",
                                              y = "Divergence")
theta.div.scatter <- theta.div.scatter + annotate("text", x = 0.0138, y = 0.085, label = "Spearman's rank partial correlation = 0.64, p = 0.0008")

ggsave("Data/iSMC/theta_paper/submission/Figures/divergence.zt.200kb.pdf", plot = theta.div.scatter, dpi = 500, width = 7, height = 7)


########################################
#
# Mutation x Evolutionary (Protein) Rates --- 50 kb
#
########################################

# loads 
zt.raw <- read.table("Data/iSMC/theta_paper/real_data/zt.misc/PerGeneAllVariables.csv",
                     sep = ",", header = T, fill = T, stringsAsFactors = T)
zt.tbl <- zt.raw[which(zt.raw$Chr == 1), c(4, 5, 6, 8, 24, 25, 29)]

zt.evol <- zt.tbl[which(zt.tbl$PiN_PiS > 0),]
zt.evol <- zt.evol[which(zt.evol$PiN_PiS < 1),]

# to get the ranges
diversity.zt <- diversity.zt.50kb[which(intersect.50kb == F),]
zt.50kb <- as.data.frame(cbind(diversity.zt[,(1:3)], zt.lands.50kb)) 
names(zt.50kb) <- c("Chr", "Start", "End", names(zt.lands.50kb))
zt.50kb$Chr <- 1







############### 
# grouping per window of constant theta
zt.lands.gr <- makeGRangesFromDataFrame(zt.50kb)
values(zt.lands.gr) <- zt.50kb[,(4:7)]
evolrate.gr <- makeGRangesFromDataFrame(zt.evol)
values(evolrate.gr) <- zt.evol[,(4:7)]

hits <- findOverlaps(zt.lands.gr, evolrate.gr) 
ranges(evolrate.gr)[subjectHits(hits)] <- ranges(zt.lands.gr)[queryHits(hits)]

evolrate.gr.df <- as.data.frame(evolrate.gr)
zt.lands.gr.df <- as.data.frame(zt.lands.gr)

# gets mean of a statistic for all genes in the window
tmp <- ddply(.data = evolrate.gr.df[-c(1,5)], .variables = "start", .fun = colMeans, .progress = "text")
evolrate.gr.df.2 <- evolrate.gr.df[!duplicated(evolrate.gr.df$start),]
evolrate.gr.df.2$score <- tmp$PiNPiS

# due to previously filtering out windows based on missing data,
# some genes do not match to a window
evolrate.gr.df.2 <- evolrate.gr.df.2[which(((evolrate.gr.df.2$width - 1) %% 50000) == 0),]

zt.lands.gr.df.2 <- zt.lands.gr.df[which(zt.lands.gr.df$start %in% evolrate.gr.df.2$start),]

zt.lands.evolrate <- merge(zt.lands.gr.df.2, evolrate.gr.df.2, by = "start")

# for filtering when score is PiN/PiS or dN/dS
zt.lands.evolrate <- zt.lands.evolrate[which(zt.lands.evolrate$score > 0),]
zt.lands.evolrate <- zt.lands.evolrate[which(zt.lands.evolrate$score < 1),]

cor.test(zt.lands.evolrate$score, zt.lands.evolrate$theta, method = "spearman")
pcor.test(zt.lands.evolrate$score, zt.lands.evolrate$theta, zt.lands.evolrate$rho, method = "spearman") 

# plot
evolrate.scatter <- ggplot(data = zt.lands.evolrate, aes(x = zt.lands.evolrate$theta, y = zt.lands.evolrate$score, color =  zt.lands.evolrate$rho))
evolrate.scatter <- evolrate.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
evolrate.scatter <- evolrate.scatter + theme(plot.title = element_text(size = 20),
                                             axis.title.x = element_text(size = 20),
                                             axis.title.y = element_text(size = 20))
evolrate.scatter <- evolrate.scatter + labs(title = NULL, # "Mutation rate vs Divergence D. melanogaster : D. yakuba",
                                            color = expression(rho),
                                            x = "Scaled Mutation Rate in D. melanogaster",
                                            y = "Evol. Rate")
evolrate.scatter




############### 
# grouping per gene coordinate
zt.lands.gr <- makeGRangesFromDataFrame(zt.50kb)
values(zt.lands.gr) <- zt.50kb[,(4:7)]
evolrate.gr <- makeGRangesFromDataFrame(zt.evol)
values(evolrate.gr) <- zt.evol[,(4:7)]

hits <- findOverlaps(evolrate.gr, zt.lands.gr) 

evolrate.gr.df <- as.data.frame(evolrate.gr[queryHits(hits)], row.names = NULL)
zt.lands.gr.df <- as.data.frame(zt.lands.gr[subjectHits(hits)], row.names = NULL)

zt.lands.evolrate <- cbind.data.frame(zt.lands.gr.df[,c(2, 3, 6:9)], evolrate.gr.df[,c(2, 3, 7:9)])
names(zt.lands.evolrate)[1] <- "start.window"
names(zt.lands.evolrate)[2] <- "end.window"
names(zt.lands.evolrate)[7] <- "start.gene"
names(zt.lands.evolrate)[8] <- "end.gene"


# basic checking
cor.test(zt.lands.evolrate$diversity, zt.lands.evolrate$PiS, method = "spearman")
cor.test(zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$PiS, method = "spearman")
cor.test(zt.lands.evolrate$diversity, zt.lands.evolrate$PiN_PiS, method = "spearman")
pcor.test(zt.lands.evolrate$diversity, zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$theta, method = "spearman")


#
cor.test(zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$theta, method = "spearman")
cor.test(zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$rho, method = "spearman")
pcor.test(zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$theta, zt.lands.evolrate$rho, method = "spearman") 
# decomposing PiN/PiS
cor.test(zt.lands.evolrate$PiS, zt.lands.evolrate$theta, method = "spearman")
cor.test(zt.lands.evolrate$PiN, zt.lands.evolrate$theta, method = "spearman")
pcor.test(zt.lands.evolrate$PiN, zt.lands.evolrate$theta, zt.lands.evolrate$rho, method = "spearman") 

# dN/dS is affected by both negative and positive selection
cor.test(zt.lands.evolrate$dS, zt.lands.evolrate$theta, method = "spearman")
cor.test(zt.lands.evolrate$dN, zt.lands.evolrate$theta, method = "spearman")
cor.test(zt.lands.evolrate$dNdS, zt.lands.evolrate$theta, method = "spearman")


# plot
evolrate.scatter <- ggplot(data = zt.lands.evolrate, aes(x = zt.lands.evolrate$theta, y = zt.lands.evolrate$dNdS, color =  zt.lands.evolrate$PiN_PiS))
evolrate.scatter <- evolrate.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
evolrate.scatter <- evolrate.scatter + theme(plot.title = element_text(size = 20),
                                             axis.title.x = element_text(size = 20),
                                             axis.title.y = element_text(size = 20))
evolrate.scatter <- evolrate.scatter + labs(title = NULL, # "Mutation rate vs Divergence D. melanogaster : D. yakuba",
                                            color = expression(pi),
                                            x = "Scaled Mutation Rate in D. melanogaster",
                                            y = "Evol. Rate")
evolrate.scatter

########################################
#
# Mutation x Evolutionary (Protein) Rates --- 200 kb
#
########################################

# loads 
zt.raw <- read.table("Data/iSMC/theta_paper/real_data/zt.misc/PerGeneAllVariables.csv",
                     sep = ",", header = T, fill = T, stringsAsFactors = T)
zt.tbl <- zt.raw[which(zt.raw$Chr == 1), c(4, 5, 6, 8, 24, 25, 29)]

zt.evol <- zt.tbl[which(zt.tbl$PiN_PiS > 0),]
zt.evol <- zt.evol[which(zt.evol$PiN_PiS < 1),]

# to get the ranges
diversity.zt <- diversity.zt.200kb[which(intersect.200kb == F),]
zt.200kb <- as.data.frame(cbind(diversity.zt[,(1:3)], zt.lands.200kb)) 
names(zt.200kb) <- c("Chr", "Start", "End", names(zt.lands.200kb))
zt.200kb$Chr <- 1







############### 
# grouping per window of constant theta
zt.lands.gr <- makeGRangesFromDataFrame(zt.200kb)
values(zt.lands.gr) <- zt.200kb[,(4:7)]
evolrate.gr <- makeGRangesFromDataFrame(zt.evol)
values(evolrate.gr) <- zt.evol[,(4:7)]

hits <- findOverlaps(zt.lands.gr, evolrate.gr) 
ranges(evolrate.gr)[subjectHits(hits)] <- ranges(zt.lands.gr)[queryHits(hits)]

evolrate.gr.df <- as.data.frame(evolrate.gr)
zt.lands.gr.df <- as.data.frame(zt.lands.gr)

# gets mean of a statistic for all genes in the window
tmp <- ddply(.data = evolrate.gr.df[-c(1,5)], .variables = "start", .fun = colMeans, .progress = "text")
evolrate.gr.df.2 <- evolrate.gr.df[!duplicated(evolrate.gr.df$start),]
evolrate.gr.df.2$score <- tmp$PiNPiS

# due to previously filtering out windows based on missing data,
# some genes do not match to a window
evolrate.gr.df.2 <- evolrate.gr.df.2[which(((evolrate.gr.df.2$width - 1) %% 200000) == 0),]

zt.lands.gr.df.2 <- zt.lands.gr.df[which(zt.lands.gr.df$start %in% evolrate.gr.df.2$start),]

zt.lands.evolrate <- merge(zt.lands.gr.df.2, evolrate.gr.df.2, by = "start")

# for filtering when score is PiN/PiS or dN/dS
zt.lands.evolrate <- zt.lands.evolrate[which(zt.lands.evolrate$score > 0),]
zt.lands.evolrate <- zt.lands.evolrate[which(zt.lands.evolrate$score < 1),]

cor.test(zt.lands.evolrate$score, zt.lands.evolrate$theta, method = "spearman")
pcor.test(zt.lands.evolrate$score, zt.lands.evolrate$theta, zt.lands.evolrate$rho, method = "spearman") 

# plot
evolrate.scatter <- ggplot(data = zt.lands.evolrate, aes(x = zt.lands.evolrate$theta, y = zt.lands.evolrate$score, color =  zt.lands.evolrate$rho))
evolrate.scatter <- evolrate.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
evolrate.scatter <- evolrate.scatter + theme(plot.title = element_text(size = 20),
                                             axis.title.x = element_text(size = 20),
                                             axis.title.y = element_text(size = 20))
evolrate.scatter <- evolrate.scatter + labs(title = NULL, # "Mutation rate vs Divergence D. melanogaster : D. yakuba",
                                            color = expression(rho),
                                            x = "Scaled Mutation Rate in D. melanogaster",
                                            y = "Evol. Rate")
evolrate.scatter




############### 
# grouping per gene coordinate
zt.lands.gr <- makeGRangesFromDataFrame(zt.200kb)
values(zt.lands.gr) <- zt.200kb[,(4:7)]
evolrate.gr <- makeGRangesFromDataFrame(zt.evol)
values(evolrate.gr) <- zt.evol[,(4:7)]

hits <- findOverlaps(evolrate.gr, zt.lands.gr) 

evolrate.gr.df <- as.data.frame(evolrate.gr[queryHits(hits)], row.names = NULL)
zt.lands.gr.df <- as.data.frame(zt.lands.gr[subjectHits(hits)], row.names = NULL)

zt.lands.evolrate <- cbind.data.frame(zt.lands.gr.df[,c(2, 3, 6:9)], evolrate.gr.df[,c(2, 3, 7:9)])
names(zt.lands.evolrate)[1] <- "start.window"
names(zt.lands.evolrate)[2] <- "end.window"
names(zt.lands.evolrate)[7] <- "start.gene"
names(zt.lands.evolrate)[8] <- "end.gene"


# basic checking
cor.test(zt.lands.evolrate$diversity, zt.lands.evolrate$PiS, method = "spearman")
cor.test(zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$PiS, method = "spearman")
cor.test(zt.lands.evolrate$diversity, zt.lands.evolrate$PiN_PiS, method = "spearman")
pcor.test(zt.lands.evolrate$diversity, zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$theta, method = "spearman")


#
cor.test(zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$theta, method = "spearman")
cor.test(zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$rho, method = "spearman")
pcor.test(zt.lands.evolrate$PiN_PiS, zt.lands.evolrate$theta, zt.lands.evolrate$rho, method = "spearman") 
# decomposing PiN/PiS
cor.test(zt.lands.evolrate$PiS, zt.lands.evolrate$theta, method = "spearman")
cor.test(zt.lands.evolrate$PiN, zt.lands.evolrate$theta, method = "spearman")
pcor.test(zt.lands.evolrate$PiN, zt.lands.evolrate$theta, zt.lands.evolrate$rho, method = "spearman") 

# dN/dS is affected by both negative and positive selection
cor.test(zt.lands.evolrate$dS, zt.lands.evolrate$theta, method = "spearman")
cor.test(zt.lands.evolrate$dN, zt.lands.evolrate$theta, method = "spearman")
cor.test(zt.lands.evolrate$dNdS, zt.lands.evolrate$theta, method = "spearman")


# plot
evolrate.scatter <- ggplot(data = zt.lands.evolrate, aes(x = zt.lands.evolrate$theta, y = zt.lands.evolrate$dNdS, color =  zt.lands.evolrate$PiN_PiS))
evolrate.scatter <- evolrate.scatter + geom_point(shape = 16, size = 4,  alpha = 0.8)
evolrate.scatter <- evolrate.scatter + theme(plot.title = element_text(size = 20),
                                             axis.title.x = element_text(size = 20),
                                             axis.title.y = element_text(size = 20))
evolrate.scatter <- evolrate.scatter + labs(title = NULL, # "Mutation rate vs Divergence D. melanogaster : D. yakuba",
                                            color = expression(pi),
                                            x = "Scaled Mutation Rate in D. melanogaster",
                                            y = "Evol. Rate")
evolrate.scatter

