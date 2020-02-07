# Created: 21/09/2018
# Last modified: 24/11/2018
# Author: Gustavo Barroso 

library(reshape2)
library(ggplot2)
library(scales) 
library(cowplot)
library(grid)
library(MASS)

setwd("~/Data/iSMC/theta_paper/simulated_data/")

num_reps <- 10
bin_sizes <- c(50000, 200e+3, 500e+3, 1e+6)

###################################################
#
# Heterogeneous Theta Landscapes in Flat Rho
#
###################################################

# simulation parameters: alpha = 0.5, g = 1e-5
## matrices to store correlations and p-values
r2.mat.sim1.alpha0.5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim1.alpha0.5_100kb) <- c("50kb.1a", "200kb.1a", "500kb.1a", "1Mb.1a")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim1.alpha0.5_100kb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_100kb <- read.table(paste("sim/flat_demography/flat_rho/maps/theta.flat_rho.alpha0.5.len100kb.",
                                         as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/flat_rho/modulated_theta/alpha_0.5/avg_length_100kb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    
    r2 <- cor.test(sim_alpha0.5_100kb$sim, tmp[,3], method = "pearson") 
    r2.mat.sim1.alpha0.5_100kb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parameters: alpha = 0.5, g = 1e-6
## matrices to store correlations and p-values
r2.mat.sim1.alpha0.5_1mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim1.alpha0.5_1mb) <- c("50kb.1b", "200kb.1b", "500kb.1b", "1Mb.1b")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim1.alpha0.5_1mb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_1mb <- read.table(paste("sim/flat_demography/flat_rho/maps/theta.flat_rho.alpha0.5.len1Mb.",
                                         as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/flat_rho/modulated_theta/alpha_0.5/avg_length_1Mb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_1mb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim1.alpha0.5_1mb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parametes: alpha = 5, g = 1e-5
## matrices to store correlations and p-values
r2.mat.sim1.alpha5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim1.alpha5_100kb) <- c("50kb.1c", "200kb.1c", "500kb.1c", "1Mb.1c")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim1.alpha5_100kb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_100kb <- read.table(paste("sim/flat_demography/flat_rho/maps/theta.flat_rho.alpha5.len100kb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/flat_rho/modulated_theta/alpha_5/avg_length_100kb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim1.alpha5_100kb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parametes alpha = 5, g = 1e-6
## matrices to store correlations and p-values
r2.mat.sim1.alpha5_1mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim1.alpha5_1mb) <- c("50kb.1d", "200kb.1d", "500kb.1d", "1Mb.1d")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim1.alpha5_1mb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_100kb <- read.table(paste("sim/flat_demography/flat_rho/maps/theta.flat_rho.alpha5.len1Mb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/flat_rho/modulated_theta/alpha_5/avg_length_1Mb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim1.alpha5_1mb[j, i] <- r2$estimate ^ 2
  }
}

# combines correlations from all combinations of parameters
r2.sim1 <- rbind(r2.mat.sim1.alpha0.5_100kb, r2.mat.sim1.alpha0.5_1mb, r2.mat.sim1.alpha5_100kb, r2.mat.sim1.alpha5_1mb)

r2.sim1$Parameters <- c(rep("A = 0.5; g = 100 kb", 4),
                         rep("A = 0.5; g = 1 Mb", 4),
                         rep("A = 5; g = 100 kb", 4),
                         rep("A = 5; g = 1 Mb", 4))
r2.sim1$bin.size <- c(rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4))
write.csv(r2.sim1, "inf/flat_demography/flat_rho/modulated_theta/r2.sim1.csv", row.names = F)


molten.cor <- melt(data = r2.sim1, id.vars = c("Parameters", "bin.size"), na.rm = T) 
sim1.model.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Parameters))
sim1.model.boxplot <- sim1.model.boxplot + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) 
sim1.model.boxplot <- sim1.model.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
sim1.model.boxplot <- sim1.model.boxplot + labs(title = NULL, x = NULL, y = expression(R^2))
sim1.model.boxplot <- sim1.model.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
sim1.model.boxplot <- sim1.model.boxplot + ylim(0.1, 1.0)
sim1.model.boxplot <- sim1.model.boxplot + theme(legend.position = "none", legend.title = NULL,
                                                 axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
#sim1.model.boxplot
ggsave(path = "inf/flat_demography/flat_rho/modulated_theta/", filename = "sim1.png", plot = sim1.model.boxplot, device = "png", dpi = 500)

###################################################
#
# Heterogeneous Rho
#
###################################################

# THETA-MODULATED (flat rho)

# simulation parameters: alpha = 0.5, g = 1e-5
## matrices to store correlations and p-values
r2.mat.sim2a.alpha0.5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim2a.alpha0.5_100kb) <- c("50kb.1a", "200kb.1a", "500kb.1a", "1Mb.1a")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim2a.alpha0.5_100kb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_100kb <- read.table(paste("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha0.5.len100kb.",
                                         as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/gamma_rho/modulated_theta/alpha_0.5/avg_length_100kb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim2a.alpha0.5_100kb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parameters: alpha = 0.5, g = 1e-6
## matrices to store correlations and p-values
r2.mat.sim2a.alpha0.5_1mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim2a.alpha0.5_1mb) <- c("50kb.1b", "200kb.1b", "500kb.1b", "1Mb.1b")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim2a.alpha0.5_1mb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_1mb <- read.table(paste("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha0.5.len1Mb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/gamma_rho/modulated_theta/alpha_0.5/avg_length_1Mb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_1mb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim2a.alpha0.5_1mb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parametes: alpha = 5, g = 1e-5
## matrices to store correlations and p-values
r2.mat.sim2a.alpha5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim2a.alpha5_100kb) <- c("50kb.1c", "200kb.1c", "500kb.1c", "1Mb.1c")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim2a.alpha5_100kb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_100kb <- read.table(paste("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha5.len100kb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/gamma_rho/modulated_theta/alpha_5/avg_length_100kb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim2a.alpha5_100kb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parametes alpha = 5, g = 1e-6
## matrices to store correlations and p-values
r2.mat.sim2a.alpha5_1mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim2a.alpha5_1mb) <- c("50kb.1d", "200kb.1d", "500kb.1d", "1Mb.1d")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim2a.alpha5_1mb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_100kb <- read.table(paste("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha5.len1Mb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/gamma_rho/modulated_theta/alpha_5/avg_length_1Mb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim2a.alpha5_1mb[j, i] <- r2$estimate ^ 2
  }
}

# combines correlations from all combinations of parameters
r2.sim2a <- rbind(r2.mat.sim2a.alpha0.5_100kb, r2.mat.sim2a.alpha0.5_1mb, r2.mat.sim2a.alpha5_100kb, r2.mat.sim2a.alpha5_1mb)

r2.sim2a$Parameters <- c(rep("A = 0.5; g = 100 kb", 4),
                         rep("A = 0.5; g = 1 Mb", 4),
                         rep("A = 5; g = 100 kb", 4),
                         rep("A = 5; g = 1 Mb", 4))
r2.sim2a$bin.size <- c(rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4))
write.csv(r2.sim2a, "inf/flat_demography/gamma_rho/modulated_theta/r2.single_mod.csv", row.names = F)

molten.cor <- melt(data = r2.sim2a, id.vars = c("Parameters", "bin.size"), na.rm = T) 
sim2a.model.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Parameters))
sim2a.model.boxplot <- sim2a.model.boxplot + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) 
sim2a.model.boxplot <- sim2a.model.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
sim2a.model.boxplot <- sim2a.model.boxplot + labs(title = NULL, x = "Bin Size", y = expression(R^2))
sim2a.model.boxplot <- sim2a.model.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
sim2a.model.boxplot <- sim2a.model.boxplot + ylim(0.1, 1.0)
sim2a.model.boxplot <- sim2a.model.boxplot + theme(legend.position = "none", #"top",
                                                 axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))



# DOUBLE-MODULATED

# simulation parameters: alpha = 0.5, g = 1e-5
## matrices to store correlations and p-values
r2.mat.sim2b.alpha0.5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim2b.alpha0.5_100kb) <- c("50kb.1a", "200kb.1a", "500kb.1a", "1Mb.1a")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim2b.alpha0.5_100kb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_100kb <- read.table(paste("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha0.5.len100kb.",
                                         as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_0.5/avg_length_100kb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim2b.alpha0.5_100kb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parameters: alpha = 0.5, g = 1e-6
## matrices to store correlations and p-values
r2.mat.sim2b.alpha0.5_1mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim2b.alpha0.5_1mb) <- c("50kb.1b", "200kb.1b", "500kb.1b", "1Mb.1b")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim2b.alpha0.5_1mb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_1mb <- read.table(paste("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha0.5.len1Mb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_0.5/avg_length_1Mb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_1mb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim2b.alpha0.5_1mb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parametes: alpha = 5, g = 1e-5
## matrices to store correlations and p-values
r2.mat.sim2b.alpha5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim2b.alpha5_100kb) <- c("50kb.1c", "200kb.1c", "500kb.1c", "1Mb.1c")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim2b.alpha5_100kb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_100kb <- read.table(paste("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha5.len100kb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_100kb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim2b.alpha5_100kb[j, i] <- r2$estimate ^ 2
  }
}

# simulation parametes alpha = 5, g = 1e-6
## matrices to store correlations and p-values
r2.mat.sim2b.alpha5_1mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim2b.alpha5_1mb) <- c("50kb.1d", "200kb.1d", "500kb.1d", "1Mb.1d")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim2b.alpha5_1mb) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_100kb <- read.table(paste("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha5.len1Mb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_1Mb/rep_", as.character(i),
                            "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim2b.alpha5_1mb[j, i] <- r2$estimate ^ 2
  }
}

# combines correlations from all combinations of parameters
r2.sim2b <- rbind(r2.mat.sim2b.alpha0.5_100kb, r2.mat.sim2b.alpha0.5_1mb, r2.mat.sim2b.alpha5_100kb, r2.mat.sim2b.alpha5_1mb)

r2.sim2b$Parameters <- c(rep("A = 0.5; g = 100 kb", 4),
                          rep("A = 0.5; g = 1 Mb", 4),
                          rep("A = 5; g = 100 kb", 4),
                          rep("A = 5; g = 1 Mb", 4))
r2.sim2b$bin.size <- c(rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4))
write.table(r2.sim2b, "inf/flat_demography/gamma_rho/multi_modulated/r2.double_mod.csv")

molten.cor <- melt(data = r2.sim2b, id.vars = c("Parameters", "bin.size"), na.rm = T) 
sim2b.model.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Parameters))
sim2b.model.boxplot <- sim2b.model.boxplot + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) 
sim2b.model.boxplot <- sim2b.model.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
sim2b.model.boxplot <- sim2b.model.boxplot + labs(title = NULL, x = "Bin Size", y = NULL)
sim2b.model.boxplot <- sim2b.model.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
sim2b.model.boxplot <- sim2b.model.boxplot + ylim(0.1, 1.0)
sim2b.model.boxplot <- sim2b.model.boxplot + theme(legend.position = "none", #"top",
                                                 axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))


# comparing accuracy of single- and double-modulated models
r2.rho <- rbind(r2.sim2a, r2.sim2b)
r2.rho$Model <- c(rep("single-modulated", 16), rep("double-modulated", 16))
write.csv(r2.rho, "inf/flat_demography/gamma_rho/r2.rho.csv", row.names = F)

molten.cor <- melt(data = r2.rho[-which(names(r2.rho) == "Parameters")], id.vars = c("Model", "bin.size")) 
rho.models.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Model))
rho.models.boxplot <- rho.models.boxplot + scale_fill_manual(values = c("#fc8d62", "#8da0cb"))
rho.models.boxplot <- rho.models.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
rho.models.boxplot <- rho.models.boxplot + labs(title = "Bin sizes", x = "Bin Size", y = "R^2")
rho.models.boxplot <- rho.models.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
rho.models.boxplot <- rho.models.boxplot + ylim(0.1, 1.0)
rho.models.boxplot <- rho.models.boxplot + theme(legend.position = "none", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
ggsave(path = "inf/flat_demography/gamma_rho/", filename = "sim2.png", plot = rho.models.boxplot, device = "png", dpi = 500)

# to facilitate putting labels in the x-axis 
# combinations of (S)mall, (L)arge, (F)requent, (R)are changes in mutation rate along the genome
r2.rho$Parameters <- rep(c(rep("L / F", 4), rep("L / R", 4), rep("S / F", 4), rep("S / R", 4)), 2)

molten.r2.2 <- melt(data = r2.rho[-which(names(r2.rho) == "bin.size")], id.vars = c("Model", "Parameters")) 
rho.models.boxplot.2 <- ggplot(data = molten.r2.2, aes(y = value, x = Parameters, fill = Model))
rho.models.boxplot.2 <- rho.models.boxplot.2 + scale_fill_manual(values = c("#fc8d62", "#8da0cb"))
rho.models.boxplot.2 <- rho.models.boxplot.2 + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
rho.models.boxplot.2 <- rho.models.boxplot.2 + labs(title = "Landscapes", x = "Change in mu", y = NULL)
rho.models.boxplot.2 <- rho.models.boxplot.2 + scale_x_discrete(limits = as.character(c("L / F", "L / R", "S / F", "S / R")))
rho.models.boxplot.2 <- rho.models.boxplot.2 + ylim(0.1, 1.0)
rho.models.boxplot.2 <- rho.models.boxplot.2 + theme(legend.position = "none", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
#rho.models.boxplot.2 <- rho.models.boxplot.2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))

legend <- get_legend(rho.models.boxplot + theme(legend.position="bottom"))

rho.fig <- plot_grid(rho.models.boxplot, rho.models.boxplot.2, nrow = 1, ncol = 2, labels = "AUTO", label_size = 18)
rho.fig <- plot_grid(rho.fig, legend, ncol = 1, rel_heights = c(1, .2))
rho.fig

cowplot::ggsave("Figure_S2.pdf", plot = rho.fig, device = "pdf", dpi = 500, width = 16, height = 8,
                path = "~/Data/iSMC/theta_paper/simulated_data/inf/flat_demography/gamma_rho/")


###################################################
#
# Maps & Scatterplots REP 1
#
###################################################


# Illustration of maps and scatterplots for different bin sizes (replicate #5), in case of jointly-inferred demography
## 50 kb
theta.50kb.bin <- read.table("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha5.len1Mb.50000.bins.txt")
rep1.50kb.bin <- read.table("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_1Mb/rep_1/theta.theta_rep_1.50kb.1-30000000.bedgraph", header = T)
theta.50kb <- as.data.frame(cbind(theta.50kb.bin$sim, rep1.50kb.bin$diploid_0, theta.50kb.bin$bin))
colnames(theta.50kb) <- c("theta.sim", "theta.inf", "bin")

theta.50kb.scatter <- ggplot(data = theta.50kb, aes(x = theta.50kb$theta.sim, y = theta.50kb$theta.inf))
theta.50kb.scatter <- theta.50kb.scatter + geom_point(shape = 16, size = 4,  alpha = 0.5)
theta.50kb.scatter <- theta.50kb.scatter + geom_smooth(method = "lm", colour = 'magenta', linetype = 2, se = F)
theta.50kb.scatter <- theta.50kb.scatter + theme(legend.position = "none",
                                             axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24))
theta.50kb.scatter <- theta.50kb.scatter + labs(title = NULL, x = "Simulation", y = "Inference")

molten.df <- melt(data = theta.50kb, id.vars = "bin")
theta.50kb.map <- ggplot(data = molten.df, aes(x = bin * 50, y = value))
theta.50kb.map <- theta.50kb.map + geom_line(data = subset(molten.df, variable == "theta.sim"), colour = "black", size = 1.2)
theta.50kb.map <- theta.50kb.map + geom_line(data = subset(molten.df, variable == "theta.inf"), colour = "salmon", size = 0.8)
theta.50kb.map <- theta.50kb.map + labs(title = NULL, x = "Position (kb)", y = expression(theta))
theta.50kb.map <- theta.50kb.map + theme(axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 36),
                                     axis.text = element_text(size = 10), legend.position = "none")
theta.50kb.map <- theta.50kb.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

## 200 kb
theta.200kb.bin <- read.table("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha5.len1Mb.2e+05.bins.txt")
rep1.200kb.bin <- read.table("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_1Mb/rep_1/theta.theta_rep_1.200kb.1-30000000.bedgraph", header = T)
theta.200kb <- as.data.frame(cbind(theta.200kb.bin$sim, rep1.200kb.bin$diploid_0, theta.200kb.bin$bin))
colnames(theta.200kb) <- c("theta.sim", "theta.inf", "bin")

theta.200kb.scatter <- ggplot(data = theta.200kb, aes(x = theta.200kb$theta.sim, y = theta.200kb$theta.inf))
theta.200kb.scatter <- theta.200kb.scatter + geom_point(shape = 16, size = 4,  alpha = 0.6)
theta.200kb.scatter <- theta.200kb.scatter + geom_smooth(method = "lm", colour = 'magenta', linetype = 2, se = F)
theta.200kb.scatter <- theta.200kb.scatter + theme(legend.position = "none",
                                               axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24))
theta.200kb.scatter <- theta.200kb.scatter + labs(title = NULL, x = "Simulation", y = " ")

molten.df <- melt(data = theta.200kb, id.vars = "bin")
theta.200kb.map <- ggplot(data = molten.df, aes(x = bin * 50, y = value))
theta.200kb.map <- theta.200kb.map + geom_line(data = subset(molten.df, variable == "theta.sim"), colour = "black", size = 1.2)
theta.200kb.map <- theta.200kb.map + geom_line(data = subset(molten.df, variable == "theta.inf"), colour = "salmon", size = 0.8)
theta.200kb.map <- theta.200kb.map + labs(title = NULL, x = "Position (kb)", y = " ")
theta.200kb.map <- theta.200kb.map + theme(axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 36),
                                       axis.text = element_text(size = 10), legend.position = "none")
theta.200kb.map <- theta.200kb.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

## 500 kb
theta.500kb.bin <- read.table("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha5.len1Mb.5e+05.bins.txt")
rep1.500kb.bin <- read.table("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_1Mb/rep_1/theta.theta_rep_1.500kb.1-30000000.bedgraph", header = T)
theta.500kb <- as.data.frame(cbind(theta.500kb.bin$sim, rep1.500kb.bin$diploid_0, theta.500kb.bin$bin))
colnames(theta.500kb) <- c("theta.sim", "theta.inf", "bin")

theta.500kb.scatter <- ggplot(data = theta.500kb, aes(x = theta.500kb$theta.sim, y = theta.500kb$theta.inf))
theta.500kb.scatter <- theta.500kb.scatter + geom_point(shape = 16, size = 4,  alpha = 0.75)
theta.500kb.scatter <- theta.500kb.scatter + geom_smooth(method = "lm", colour = 'magenta', linetype = 2, se = F)
theta.500kb.scatter <- theta.500kb.scatter + theme(legend.position = "none",
                                               axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24))
theta.500kb.scatter <- theta.500kb.scatter + labs(title = NULL, x = "Simulation", y = " ")

molten.df <- melt(data = theta.500kb, id.vars = "bin")
theta.500kb.map <- ggplot(data = molten.df, aes(x = bin * 50, y = value))
theta.500kb.map <- theta.500kb.map + geom_line(data = subset(molten.df, variable == "theta.sim"), colour = "black", size = 1.2)
theta.500kb.map <- theta.500kb.map + geom_line(data = subset(molten.df, variable == "theta.inf"), colour = "salmon", size = 0.8)
theta.500kb.map <- theta.500kb.map + labs(title = NULL, x = "Position (kb)", y = " ")
theta.500kb.map <- theta.500kb.map + theme(axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 36),
                                       axis.text = element_text(size = 10), legend.position = "none")
theta.500kb.map <- theta.500kb.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

## 1 Mb
theta.1mb.bin <- read.table("sim/flat_demography/gamma_rho/maps/theta.gamma_rho.alpha5.len1Mb.1e+06.bins.txt")
rep1.1mb.bin <- read.table("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_1Mb/rep_1/theta.theta_rep_1.1mb.1-30000000.bedgraph", header = T)
theta.1mb <- as.data.frame(cbind(theta.1mb.bin$sim, rep1.1mb.bin$diploid_0, theta.1mb.bin$bin))
colnames(theta.1mb) <- c("theta.sim", "theta.inf", "bin")

theta.1mb.scatter <- ggplot(data = theta.1mb, aes(x = theta.1mb$theta.sim, y = theta.1mb$theta.inf))
theta.1mb.scatter <- theta.1mb.scatter + geom_point(shape = 16, size = 4,  alpha = 0.9)
theta.1mb.scatter <- theta.1mb.scatter + geom_smooth(method = "lm", colour = 'magenta', linetype = 2, se = F)
theta.1mb.scatter <- theta.1mb.scatter + theme(axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22))
theta.1mb.scatter <- theta.1mb.scatter + labs(title = NULL, x = "Simulation", y = "Inference")

molten.df <- melt(data = theta.1mb, id.vars = "bin")
theta.1mb.map <- ggplot(data = molten.df, aes(x = bin * 50, y = value))
theta.1mb.map <- theta.1mb.map + geom_line(data = subset(molten.df, variable == "theta.sim"), colour = "black", size = 1.2)
theta.1mb.map <- theta.1mb.map + geom_line(data = subset(molten.df, variable == "theta.inf"), colour = "salmon", size = 0.8)
theta.1mb.map <- theta.1mb.map + labs(title = NULL, x = "Position (kb)", y = expression(theta))
theta.1mb.map <- theta.1mb.map + theme(axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 36))
theta.1mb.map <- theta.1mb.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

Fig.3 <- plot_grid(theta.50kb.scatter, theta.200kb.scatter, theta.500kb.scatter, theta.1mb.scatter,
                   theta.50kb.map, theta.200kb.map, theta.500kb.map, theta.1mb.map,
                   nrow = 2, ncol = 4, labels = "AUTO", label_size = 36, scale = 0.9)
Fig.3
setTimeLimit(cpu = Inf, elapsed = Inf, transient = T)
cowplot::ggsave("../submission/Figures/Figure5.pdf", plot = Fig.3, device = "pdf", dpi = 500, width = 24, height = 12)

defence_2 <- plot_grid(theta.1mb.scatter, theta.1mb.map, nrow = 1, ncol = 2, labels = NULL, scale = 0.9)

ggsave("../submission/Figures/Defence2.pdf", plot = defence_2, device = "pdf", dpi = 500, width = 16, height = 8)


###################################################
#
# Population Size Changes
#
###################################################

# matrices for correlations and p-values when demography is mis-specified as flat ('fit')
r2.mat.sim3.fix <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim3.fix) <- c("50kb.fix", "200kb.fix", "500kb.fix", "1Mb.fix")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim3.fix) <- reps

# matrices for correlations and p-values when demography is jointly inferred with the mutation map ('fit')
r2.mat.sim3.fit <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.sim3.fit) <- c("50kb.fit", "200kb.fit", "500kb.fit", "1Mb.fit")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.sim3.fit) <- reps

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters

# simulation parameters: alpha = 0.5, g = 1e-5
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_100kb <- read.table(paste("sim/bottleneck_0.5/flat_rho/maps/theta.bottleneck.alpha0.5.len100kb.",
                                         as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/bottleneck/flat_rho/fix_demo/alpha_0.5/avg_length_100kb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim3.fix[j, i] <- r2$estimate ^ 2
  }
}
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_100kb <- read.table(paste("sim/bottleneck_0.5/flat_rho/maps/theta.bottleneck.alpha0.5.len100kb.",
                                         as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/bottleneck/flat_rho/fit_demo/alpha_0.5/avg_length_100kb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim3.fit[j, i] <- r2$estimate ^ 2
  }
}
# combines matrices for the recent expansion scenarion
r2.demo.alpha0.5len100kb <- rbind(r2.mat.sim3.fix, r2.mat.sim3.fit)

# simulation parameters: alpha = 0.5, g = 1e-6
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_1Mb <- read.table(paste("sim/bottleneck_0.5/flat_rho/maps/theta.bottleneck.alpha0.5.len1Mb.",
                                         as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/bottleneck/flat_rho/fix_demo/alpha_0.5/avg_length_1Mb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_1Mb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim3.fix[j, i] <- r2$estimate ^ 2
  }
}
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha0.5_1Mb <- read.table(paste("sim/bottleneck_0.5/flat_rho/maps/theta.bottleneck.alpha0.5.len1Mb.",
                                         as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/bottleneck/flat_rho/fit_demo/alpha_0.5/avg_length_1Mb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha0.5_1Mb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim3.fit[j, i] <- r2$estimate ^ 2
  }
}
# combines matrices for the recent expansion scenarios
r2.demo.alpha0.5len1Mb <- rbind(r2.mat.sim3.fix, r2.mat.sim3.fit)

# simulation parameters: alpha = 0.5, g = 1e-5
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_100kb <- read.table(paste("sim/bottleneck_0.5/flat_rho/maps/theta.bottleneck.alpha5.len100kb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/bottleneck/flat_rho/fix_demo/alpha_5/avg_length_100kb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim3.fix[j, i] <- r2$estimate ^ 2
  }
}
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_100kb <- read.table(paste("sim/bottleneck_0.5/flat_rho/maps/theta.bottleneck.alpha5.len100kb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/bottleneck/flat_rho/fit_demo/alpha_5/avg_length_100kb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_100kb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim3.fit[j, i] <- r2$estimate ^ 2
  }
}
# combines matrices for the recent expansion scenarios
r2.demo.alpha5len100kb <- rbind(r2.mat.sim3.fix, r2.mat.sim3.fit)

# simulation parameters: alpha = 5, g = 1e-6
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_1Mb <- read.table(paste("sim/bottleneck_0.5/flat_rho/maps/theta.bottleneck.alpha5.len1Mb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/bottleneck/flat_rho/fix_demo/alpha_5/avg_length_1Mb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_1Mb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim3.fix[j, i] <- r2$estimate ^ 2
  }
}
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_alpha5_1Mb <- read.table(paste("sim/bottleneck_0.5/flat_rho/maps/theta.bottleneck.alpha5.len1Mb.",
                                       as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replicates
    tmp <- read.table(paste("inf/bottleneck/flat_rho/fit_demo/alpha_5/avg_length_1Mb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_alpha5_1Mb$sim, tmp[,3], method = "pearson")  
    r2.mat.sim3.fit[j, i] <- r2$estimate ^ 2
  }
}
# combines matrices for the recent expansion scenarios
r2.demo.alpha5len1Mb <- rbind(r2.mat.sim3.fix, r2.mat.sim3.fit)

# combines matrices
r2.sim3 <- rbind(r2.demo.alpha0.5len100kb, r2.demo.alpha0.5len1Mb, r2.demo.alpha5len100kb, r2.demo.alpha5len1Mb)
r2.sim3$Demography <- rep(c(rep("Mis-specified", 4), rep("Inferred", 4)), 4)
r2.sim3$bin.size <- rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4)

# plots all landscapes combined (comparing fix vs fit models)
molten.cor <- melt(data = r2.sim3, id.vars = c("Demography", "bin.size")) 
sim3.model.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Demography))
#sim3.model.boxplot <- sim3.model.boxplot + scale_fill_manual(values = c("#fc8d62", "#8da0cb"))
sim3.model.boxplot <- sim3.model.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
sim3.model.boxplot <- sim3.model.boxplot + labs(title = NULL, x = "Bin Size", y = expression(R^2))
sim3.model.boxplot <- sim3.model.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
sim3.model.boxplot <- sim3.model.boxplot + ylim(0.1, 1.0)
sim3.model.boxplot <- sim3.model.boxplot + theme(legend.position = "none", #"top",
                                                 axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
ggsave("theta.demography.pdf", plot = sim3.model.boxplot, device = "pdf")

r2.sim3$Parameters <- c(rep("A = 0.5, g = 1e-5", 8), rep("A = 0.5, g = 1e-6", 8),
                         rep("A = 5, g = 1e-5", 8), rep("A = 5, g = 1e-6", 8))
write.csv(r2.sim3, "inf/bottleneck/flat_rho/r2.demography.csv", row.names = F, col.names = T)

# p.values
pval.sim3 <- rbind(pval.demo.alpha0.5len100kb, pval.demo.alpha0.5len1Mb, pval.demo.alpha5len100kb, pval.demo.alpha5len1Mb)
pval.sim3$Demography <- rep(c(rep("Mis-specified", 4), rep("Inferred", 4)), 4)
pval.sim3$bin.size <- rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4)
r2.sim3$Parameters <- c(rep("A = 0.5, g = 1e-5", 8), rep("A = 0.5, g = 1e-6", 8),
                         rep("A = 5, g = 1e-5", 8), rep("A = 5, g = 1e-6", 8))
write.csv(pval.sim3, "inf/bottleneck/flat_rho/pvalues.csv", row.names = F, col.names = T)

# compares fix and fix models of demography for particular landscapes
molten.cor <- melt(data = r2.sim3[-which(names(r2.sim3) == "Parameters")], id.vars = c("Demography", "bin.size")) 
sim3.demo.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Demography))
sim3.demo.boxplot <- sim3.demo.boxplot + scale_fill_manual(values = c("#fc8d62", "#8da0cb"))
sim3.demo.boxplot <- sim3.demo.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                      binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                      alpha = 0.5)
sim3.demo.boxplot <- sim3.demo.boxplot + labs(title = NULL, x = NULL, y = expression(R^2))
sim3.demo.boxplot <- sim3.demo.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
sim3.demo.boxplot <- sim3.demo.boxplot + ylim(0.1, 1.0)
sim3.demo.boxplot <- sim3.demo.boxplot + theme(legend.position = "none", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
ggsave(path = "inf/bottleneck/flat_rho/", filename = "demo.png", plot = sim3.demo.boxplot, device = "png", dpi = 500)

# mis-specified demography
r2.sim3.mis <- r2.sim3[which(r2.sim3$Demography == "Mis-specified"),]
r2.sim3.mis <- r2.sim3.mis[-which(names(r2.sim3.mis) == "Demography")]
molten.cor <- melt(data = r2.sim3.mis, id.vars = c("Parameters", "bin.size")) 
sim3.fix.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Parameters))
sim3.fix.boxplot <- sim3.fix.boxplot + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"))
sim3.fix.boxplot <- sim3.fix.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
sim3.fix.boxplot <- sim3.fix.boxplot + labs(title = NULL, x = NULL, y = expression(R^2))
sim3.fix.boxplot <- sim3.fix.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
sim3.fix.boxplot <- sim3.fix.boxplot + ylim(0.1, 1.0)
sim3.fix.boxplot <- sim3.fix.boxplot + theme(legend.position = "none", #"top",
                                             axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))


# inferred demography
r2.sim3.inf <- r2.sim3[which(r2.sim3$Demography == "Inferred"),]
r2.sim3.inf <- r2.sim3.inf[-which(names(r2.sim3.inf) == "Demography")]
molten.cor <- melt(data = r2.sim3.inf, id.vars = c("Parameters", "bin.size")) 
sim3.fit.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Parameters))
sim3.fit.boxplot <- sim3.fit.boxplot + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"))
sim3.fit.boxplot <- sim3.fit.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                    binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                    alpha = 0.5)
sim3.fit.boxplot <- sim3.fit.boxplot + labs(title = NULL, x = NULL, y = NULL)
sim3.fit.boxplot <- sim3.fit.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
sim3.fit.boxplot <- sim3.fit.boxplot + ylim(0.1, 1.0)
sim3.fit.boxplot <- sim3.fit.boxplot + theme(legend.position = "none", #"top",
                                             axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))


###################################################
#
# Introgression
#
###################################################

# matrices for correlations and p-values when introgression happened 0.125 coal. time units ago
r2.mat.admix_0.125 <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.admix_0.125) <- c("50kb.0125", "200kb.0125", "500kb.0125", "1Mb.0125")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.admix_0.125) <- reps

# simulation parameters: alpha = 0.5, g = 1e-5
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_admix <- read.table(paste("sim/introgression/admix_time_0.125/alpha_0.5/avg_length_100kb/theta.admix.", as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replciates
    tmp <- read.table(paste("inf/introgression/admix_time_0.125/alpha_0.5/avg_length_100kb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_admix$sim, tmp[,3], method = "pearson")  
    r2.mat.admix_0.125[j, i] <- r2$estimate ^ 2
  }
}
r2.admix.0.125.alpha0.5len100kb <- r2.mat.admix_0.125

# simulation parameters: alpha = 0.5, g = 1e-6
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_admix <- read.table(paste("sim/introgression/admix_time_0.125/alpha_0.5/avg_length_1Mb/theta.admix.", as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replciates
    tmp <- read.table(paste("inf/introgression/admix_time_0.125/alpha_0.5/avg_length_1Mb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_admix$sim, tmp[,3], method = "pearson")  
    r2.mat.admix_0.125[j, i] <- r2$estimate ^ 2
  }
}
r2.admix.0.125.alpha0.5len1Mb <- r2.mat.admix_0.125

# simulation parameters: alpha = 5, g = 1e-5
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_admix <- read.table(paste("sim/introgression/admix_time_0.125/alpha_5/avg_length_100kb/theta.admix.", as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replciates
    tmp <- read.table(paste("inf/introgression/admix_time_0.125/alpha_5/avg_length_100kb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_admix$sim, tmp[,3], method = "pearson")  
    r2.mat.admix_0.125[j, i] <- r2$estimate ^ 2
  }
}
r2.admix.0.125.alpha5len100kb <- r2.mat.admix_0.125

# simulation parameters: alpha = 5, g = 1e-6
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_admix <- read.table(paste("sim/introgression/admix_time_0.125/alpha_5/avg_length_1Mb/theta.admix.", as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replciates
    tmp <- read.table(paste("inf/introgression/admix_time_0.125/alpha_5/avg_length_1Mb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_admix$sim, tmp[,3], method = "pearson")  
    r2.mat.admix_0.125[j, i] <- r2$estimate ^ 2
  }
}
r2.admix.0.125.alpha5len1Mb <- r2.mat.admix_0.125

r2.admix.0.125 <- rbind(r2.admix.0.125.alpha0.5len100kb, r2.admix.0.125.alpha0.5len1Mb,
                         r2.admix.0.125.alpha5len100kb, r2.admix.0.125.alpha5len1Mb)

r2.admix.0.125$Admix.Time <- rep("0.125 CTA", 16)
r2.admix.0.125$bin.size <- rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4)
r2.admix.0.125$Parameters <- c(rep("A = 0.5; g = 100 kb", 4),
                                rep("A = 0.5; g = 1 Mb", 4),
                                rep("A = 5; g = 100 kb", 4),
                                rep("A = 5; g = 1 Mb", 4))
write.csv(r2.admix.0.125, "inf/introgression/admix_time_0.125/r2.admix.0.125.csv", sep = "\t", row.names = F, col.names = T)

molten.cor <- melt(data = r2.admix.0.125[-which(names(r2.admix.0.125) == "Admix.Time")], id.vars = c("Parameters", "bin.size")) 
admix.0.125.model.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Parameters))
admix.0.125.model.boxplot <- admix.0.125.model.boxplot + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) 
admix.0.125.model.boxplot <- admix.0.125.model.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
admix.0.125.model.boxplot <- admix.0.125.model.boxplot + labs(title = NULL, x = NULL, y = expression(R^2))
admix.0.125.model.boxplot <- admix.0.125.model.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
admix.0.125.model.boxplot <- admix.0.125.model.boxplot + ylim(0.1, 1.0)
admix.0.125.model.boxplot <- admix.0.125.model.boxplot + theme(legend.position = "none", 
                                                               axis.title.x = element_text(size = 14),
                                                               axis.title.y = element_text(size = 14))
ggsave("inf/introgression/admix_time_0.125/admix.0.125.pdf", plot = admix.0.125.model.boxplot, device = "pdf")




# matrices for correlations and p-values when introgression happened 0.25 coal. time units ago
r2.mat.admix_0.25 <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(r2.mat.admix_0.25) <- c("50kb.025", "200kb.025", "500kb.025", "1Mb.025")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(r2.mat.admix_0.25) <- reps

# simulation parameters: alpha = 0.5, g = 1e-5
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_admix <- read.table(paste("sim/introgression/admix_time_0.25/alpha_0.5/avg_length_100kb/theta.admix.", as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replciates
    tmp <- read.table(paste("inf/introgression/admix_time_0.25/alpha_0.5/avg_length_100kb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_admix$sim, tmp[,3], method = "pearson")  
    r2.mat.admix_0.25[j, i] <- r2$estimate ^ 2
  }
}
r2.admix.0.25.alpha0.5len100kb <- r2.mat.admix_0.25

# simulation parameters: alpha = 0.5, g = 1e-6
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_admix <- read.table(paste("sim/introgression/admix_time_0.25/alpha_0.5/avg_length_1Mb/theta.admix.", as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replciates
    tmp <- read.table(paste("inf/introgression/admix_time_0.25/alpha_0.5/avg_length_1Mb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_admix$sim, tmp[,3], method = "pearson")  
    r2.mat.admix_0.25[j, i] <- r2$estimate ^ 2
  }
}
r2.admix.0.25.alpha0.5len1Mb <- r2.mat.admix_0.25

# simulation parameters: alpha = 5, g = 1e-5
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_admix <- read.table(paste("sim/introgression/admix_time_0.25/alpha_5/avg_length_100kb/theta.admix.", as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replciates
    tmp <- read.table(paste("inf/introgression/admix_time_0.25/alpha_5/avg_length_100kb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_admix$sim, tmp[,3], method = "pearson")  
    r2.mat.admix_0.25[j, i] <- r2$estimate ^ 2
  }
}
r2.admix.0.25.alpha5len100kb <- r2.mat.admix_0.25

# simulation parameters: alpha = 5, g = 1e-6
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  sim_admix <- read.table(paste("sim/introgression/admix_time_0.25/alpha_5/avg_length_1Mb/theta.admix.", as.character(bin_sizes[j]), ".bins.txt", sep = ""))
  for(i in 1:10) { # 10 replciates
    tmp <- read.table(paste("inf/introgression/admix_time_0.25/alpha_5/avg_length_1Mb/rep_", as.character(i), "/theta.theta_rep_",
                            as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    r2 <- cor.test(sim_admix$sim, tmp[,3], method = "pearson")  
    r2.mat.admix_0.25[j, i] <- r2$estimate ^ 2
  }
}
r2.admix.0.25.alpha5len1Mb <- r2.mat.admix_0.25

r2.admix.0.25 <- rbind(r2.admix.0.25.alpha0.5len100kb, r2.admix.0.25.alpha0.5len1Mb,
                        r2.admix.0.25.alpha5len100kb, r2.admix.0.25.alpha5len1Mb)

r2.admix.0.25$Admix.Time <- rep("0.25 CTA", 16)
r2.admix.0.25$bin.size <- rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4)
r2.admix.0.25$Parameters <- c(rep("A = 0.5; g = 100 kb", 4),
                               rep("A = 0.5; g = 1 Mb", 4),
                               rep("A = 5; g = 100 kb", 4),
                               rep("A = 5; g = 1 Mb", 4))
write.csv(r2.admix.0.25, "inf/introgression/admix_time_0.25/r2.admix.0.25.csv", sep = "\t", row.names = F, col.names = T)

molten.cor <- melt(data = r2.admix.0.25[-which(names(r2.admix.0.25) == "Admix.Time")], id.vars = c("Parameters", "bin.size")) 
admix.0.25.model.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Parameters))
admix.0.25.model.boxplot <- admix.0.25.model.boxplot + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) 
admix.0.25.model.boxplot <- admix.0.25.model.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                        binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                        alpha = 0.5)
admix.0.25.model.boxplot <- admix.0.25.model.boxplot + labs(title = NULL, x = NULL, y = NULL)
admix.0.25.model.boxplot <- admix.0.25.model.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
admix.0.25.model.boxplot <- admix.0.25.model.boxplot + ylim(0.1, 1.0)
admix.0.25.model.boxplot <- admix.0.25.model.boxplot + theme(legend.position = "none", 
                                                             axis.title.x = element_text(size = 14),
                                                             axis.title.y = element_text(size = 14))
ggsave("inf/introgression/admix_time_0.25/admix.0.25.pdf", plot = admix.0.25.model.boxplot, device = "pdf")


# comparing distributions for different times since introgression (in CTA)
r2.intro <- rbind(r2.admix.0.125, r2.admix.0.25)
write.csv(r2.intro, "inf/introgression/r2.introgression.csv", row.names = F)

molten.cor <- melt(data = r2.intro[-which(names(r2.intro) == "Parameters")], id.vars = c("Admix.Time", "bin.size")) 
admix.time.model.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Admix.Time))
admix.time.model.boxplot <- admix.time.model.boxplot + scale_fill_manual(values = c("#fc8d62", "#8da0cb"))
admix.time.model.boxplot <- admix.time.model.boxplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                                    binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                                    alpha = 0.5)
admix.time.model.boxplot <- admix.time.model.boxplot + labs(title = NULL, x = "Bin Size", y = NULL)
admix.time.model.boxplot <- admix.time.model.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
admix.time.model.boxplot <- admix.time.model.boxplot + ylim(0.1, 1.0)
admix.time.model.boxplot <- admix.time.model.boxplot + theme(legend.position = "top", 
                                                             axis.title.x = element_text(size = 14),
                                                             axis.title.y = element_text(size = 14))
ggsave("inf/introgression/compare.CTA.pdf", plot = admix.time.model.boxplot, device = "pdf")

molten.r2.2 <- melt(data = r2.intro[-c(12, 13)], id.vars = "Admix.Time") 
admix.time.model.boxplot.2 <- ggplot(data = molten.r2.2, aes(y = value, x = Admix.Time))
admix.time.model.boxplot.2 <- admix.time.model.boxplot.2 + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                                                    binwidth = 0.005, dotsize = 10, position = position_dodge(width=1.0),
                                                                    alpha = 0.5)
admix.time.model.boxplot.2 <- admix.time.model.boxplot.2 + labs(title = NULL, x = "Bin Size", y = NULL)
admix.time.model.boxplot.2 <- admix.time.model.boxplot.2 + scale_x_discrete(limits = as.character(c("0.125 CTA", "0.25 CTA")))
admix.time.model.boxplot.2 <- admix.time.model.boxplot.2 + ylim(0.1, 1.0)
admix.time.model.boxplot.2 <- admix.time.model.boxplot.2 + theme(legend.position = "top", 
                                                             axis.title.x = element_text(size = 14),
                                                             axis.title.y = element_text(size = 14))
ggsave("inf/introgression/compare.CTA.2.pdf", plot = admix.time.model.boxplot.2, device = "pdf")


###################################################
#
# Combining plots
#
###################################################

p <- plot_grid(sim1.model.boxplot, NULL,
               sim3.fix.boxplot, sim3.fit.boxplot,
               admix.0.125.model.boxplot, admix.0.25.model.boxplot,
               sim2a.model.boxplot, sim2b.model.boxplot,
               nrow = 4, ncol = 2, labels = c('A', ' ', 'B', 'C', 'D', 'E', 'F', 'G'), label_size = 18, scale = 0.9)

legend <- get_legend(sim1.model.boxplot + theme(legend.position="bottom"))
theme_set(theme_cowplot(font_size=12))

Fig.3 <- plot_grid(p, legend, ncol = 1, rel_widths = c(1, 0.2), rel_heights = c(1, 0.2))

Fig.3

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)
cowplot::ggsave("Figure3.pdf", plot = Fig.3, device = "pdf", dpi = 500, width = 12, height = 15,
                path = "~/Data/iSMC/theta_paper/submission/Figures/")


###################################################
#
# Correlation between theta and rho landscape in simulated data
#
###################################################


# simulation parameters: alpha = 0.5, g = 1e-5
## matrices to store correlations and p-values
cor.mat.rho.theta.alpha0.5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(cor.mat.rho.theta.alpha0.5_100kb) <- c("50kb.1a", "200kb.1a", "500kb.1a", "1Mb.1a")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(cor.mat.rho.theta.alpha0.5_100kb) <- reps
pval.mat.rho.theta.alpha0.5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(pval.mat.rho.theta.alpha0.5_100kb) <- rownames(cor.mat.rho.theta.alpha0.5_100kb)
colnames(pval.mat.rho.theta.alpha0.5_100kb) <- colnames(cor.mat.rho.theta.alpha0.5_100kb)

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  for(i in 1:10) { # 10 replicates
    tmp.theta <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_0.5/avg_length_100kb/rep_", as.character(i),
                                  "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    tmp.rho <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_0.5/avg_length_100kb/rep_", as.character(i),
                                "/rho.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    
    r2 <- cor.test(tmp.rho[,3], tmp.theta[,3], method = "pearson")  
    cor.mat.rho.theta.alpha0.5_100kb[j, i] <- r2$estimate
    pval.mat.rho.theta.alpha0.5_100kb[j, i] <- r2$p.val
  }
}

# simulation parameters: alpha = 0.5, g = 1e-6
## matrices to store correlations and p-values
cor.mat.rho.theta.alpha0.5_1Mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(cor.mat.rho.theta.alpha0.5_1Mb) <- c("50kb.1b", "200kb.1b", "500kb.1b", "1Mb.1b")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(cor.mat.rho.theta.alpha0.5_1Mb) <- reps
pval.mat.rho.theta.alpha0.5_1Mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(pval.mat.rho.theta.alpha0.5_1Mb) <- rownames(cor.mat.rho.theta.alpha0.5_1Mb)
colnames(pval.mat.rho.theta.alpha0.5_1Mb) <- colnames(cor.mat.rho.theta.alpha0.5_1Mb)

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  for(i in 1:10) { # 10 replicates
    tmp.theta <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_0.5/avg_length_1Mb/rep_", as.character(i),
                                  "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    tmp.rho <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_0.5/avg_length_1Mb/rep_", as.character(i),
                                "/rho.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    
    r2 <- cor.test(tmp.rho[,3], tmp.theta[,3], method = "pearson")  
    cor.mat.rho.theta.alpha0.5_1Mb[j, i] <- r2$estimate
    pval.mat.rho.theta.alpha0.5_1Mb[j, i] <- r2$p.val
  }
}

# simulation parameters: alpha = 5, g = 1e-5
## matrices to store correlations and p-values
cor.mat.rho.theta.alpha5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(cor.mat.rho.theta.alpha5_100kb) <- c("50kb.1c", "200kb.1c", "500kb.1c", "100kb.1c")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(cor.mat.rho.theta.alpha5_100kb) <- reps
pval.mat.rho.theta.alpha5_100kb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(pval.mat.rho.theta.alpha5_100kb) <- rownames(cor.mat.rho.theta.alpha5_100kb)
colnames(pval.mat.rho.theta.alpha5_100kb) <- colnames(cor.mat.rho.theta.alpha5_100kb)

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  for(i in 1:10) { # 10 replicates
    tmp.theta <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_100kb/rep_", as.character(i),
                                  "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    tmp.rho <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_100kb/rep_", as.character(i),
                                "/rho.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    
    r2 <- cor.test(tmp.rho[,3], tmp.theta[,3], method = "pearson")  
    cor.mat.rho.theta.alpha5_100kb[j, i] <- r2$estimate
    pval.mat.rho.theta.alpha5_100kb[j, i] <- r2$p.val
  }
}

# simulation parameters: alpha = 5, g = 1e-6
## matrices to store correlations and p-values
cor.mat.rho.theta.alpha5_1Mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(cor.mat.rho.theta.alpha5_1Mb) <- c("50kb.1d", "200kb.1d", "500kb.1d", "1Mb.1d")
reps <- character(length = num_reps)
for(i in 1: num_reps) {
  reps[i] <- paste("rep", as.character(i), sep = ".")
}
colnames(cor.mat.rho.theta.alpha5_1Mb) <- reps
pval.mat.rho.theta.alpha5_1Mb <- as.data.frame(matrix(ncol = 10, nrow = 4))
rownames(pval.mat.rho.theta.alpha5_1Mb) <- rownames(cor.mat.rho.theta.alpha5_1Mb)
colnames(pval.mat.rho.theta.alpha5_1Mb) <- colnames(cor.mat.rho.theta.alpha5_1Mb)

# these loops exploit an organised directory tree to load the correct maps for a given combination of parameters
for(j in 1:length(bin_sizes)) {
  bin_str <- paste(as.character(bin_sizes[j] / 1e+3), "kb", sep = "")
  if(j == 4) { bin_str <- "1mb" }
  for(i in 1:10) { # 10 replicates
    tmp.theta <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_1Mb/rep_", as.character(i),
                                  "/theta.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    tmp.rho <- read.table(paste("inf/flat_demography/gamma_rho/multi_modulated/alpha_5/avg_length_1Mb/rep_", as.character(i),
                                "/rho.theta_rep_", as.character(i), ".", bin_str, ".1-30000000.bedgraph", sep = ""), header = T)
    
    r2 <- cor.test(tmp.rho[,3], tmp.theta[,3], method = "pearson")  
    cor.mat.rho.theta.alpha5_1Mb[j, i] <- r2$estimate
    pval.mat.rho.theta.alpha5_1Mb[j, i] <- r2$p.val
  }
}


# combines correlations from all combinations of parameters
cor.rt <- rbind(cor.mat.rho.theta.alpha0.5_100kb, cor.mat.rho.theta.alpha0.5_1Mb, cor.mat.rho.theta.alpha5_100kb, cor.mat.rho.theta.alpha5_1Mb)
pval.rt <- rbind(pval.mat.rho.theta.alpha0.5_100kb, pval.mat.rho.theta.alpha0.5_1Mb, pval.mat.rho.theta.alpha5_100kb, pval.mat.rho.theta.alpha5_1Mb)

cor.rt$Parameters <- c(rep("A = 0.5; g = 100 kb", 4),
                       rep("A = 0.5; g = 1 Mb", 4),
                       rep("A = 5; g = 100 kb", 4),
                       rep("A = 5; g = 1 Mb", 4))
cor.rt$bin.size <- c(rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4))
write.csv(cor.rt, "inf/flat_demography/gamma_rho/multi_modulated/cor.rho.theta.csv", row.names = F)

log.pval <- log10(pval.rt)
log.pval$Parameters <- c(rep("A = 0.5; g = 100 kb", 4),
                         rep("A = 0.5; g = 1 Mb", 4),
                         rep("A = 5; g = 100 kb", 4),
                         rep("A = 5; g = 1 Mb", 4))
log.pval$bin.size <- c(rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4))

pval.rt$Parameters <- c(rep("A = 0.5; g = 100 kb", 4),
                          rep("A = 0.5; g = 1 Mb", 4),
                          rep("A = 5; g = 100 kb", 4),
                          rep("A = 5; g = 1 Mb", 4))
pval.rt$bin.size <- c(rep(c("50 kb", "200 kb", "500 kb", "1 Mb"), 4))
write.csv(pval.rt, "inf/flat_demography/gamma_rho/multi_modulated/pval.rho.theta.csv", row.names = F)



molten.cor <- melt(data = cor.rt, id.vars = c("Parameters", "bin.size"), na.rm = T) 
rho.theta.boxplot <- ggplot(data = molten.cor, aes(y = value, x = bin.size, fill = Parameters))
rho.theta.boxplot <- rho.theta.boxplot + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) 
rho.theta.boxplot <- rho.theta.boxplot + geom_boxplot(outlier.colour = "black", outlier.shape = 16, outlier.size = 1, notch=FALSE)
rho.theta.boxplot <- rho.theta.boxplot + labs(title = "Theta vs Rho (Sim.)", x = "Bin Size", y = "R^2")
rho.theta.boxplot <- rho.theta.boxplot + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
rho.theta.boxplot <- rho.theta.boxplot + theme(legend.position = "none", legend.title = NULL,
                                               axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))


molten.pval <- melt(data = log.pval, id.vars = c("Parameters", "bin.size"), na.rm = T) 
rho.theta.boxplot.2 <- ggplot(data = molten.pval, aes(y = value, x = bin.size, fill = Parameters))
rho.theta.boxplot.2 <- rho.theta.boxplot.2 + scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) 
rho.theta.boxplot.2 <- rho.theta.boxplot.2 + geom_boxplot(outlier.colour = "black", outlier.shape = 16, outlier.size = 1, notch=FALSE)
rho.theta.boxplot.2 <- rho.theta.boxplot.2 + geom_hline(yintercept = -3, colour = "magenta", linetype = "dashed")
rho.theta.boxplot.2 <- rho.theta.boxplot.2 + labs(title = "Theta vs Rho (Sim.)", x = "Bin Size", y = "Log(P-value)")
rho.theta.boxplot.2 <- rho.theta.boxplot.2 + scale_x_discrete(limits = as.character(c("1 Mb", "500 kb", "200 kb", "50 kb")))
rho.theta.boxplot.2 <- rho.theta.boxplot.2 + theme(legend.position = "none", legend.title = NULL,
                                                   axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))


p <- plot_grid(rho.theta.boxplot, rho.theta.boxplot.2,
               nrow = 1, ncol = 2, labels = "AUTO", label_size = 18)

theme_set(theme_cowplot(font_size=12))
legend <- get_legend(rho.theta.boxplot + theme(legend.position="bottom"))
Fig.S3 <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, .2))
Fig.S3

setTimeLimit(cpu = Inf, elapsed = Inf, transient = F)
cowplot::ggsave("Figure_S3.pdf", plot = Fig.S3, device = "pdf", dpi = 500, width = 12, height = 10,
                path = "~/Data/iSMC/theta_paper/simulated_data/inf/flat_demography/gamma_rho/multi_modulated/")
