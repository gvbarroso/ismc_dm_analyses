###### Li and Durbin time interval


setwd("~")
setwd("../../data2/gvbarroso/Data/iSMC/rigourous_sim/results/tmax/")


###################################################
#
# Param estimates
#
###################################################

# parameters used in the simulation
sim_params <- read.table("../../raw_data/Misc/sim_params.txt")
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
no.legend <- theme(legend.position = "none", legend.title = NULL)


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
sim.rho.50k <- read.table("../../raw_data/rig.sims.rho.50000.bins.txt")

nbins <- 3e+7 / 5e+4

# 45 single pairs of genomes
r2.vec <- numeric(length = 45)
rho.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".rho.50kb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.rho.50k$sim, tmp$diploid_1, method = "pearson")$estimate ^ 2
  rho.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.50k <- as.data.frame(r2.vec)
r2.50k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.rho.50kb.bedgraph", sep = ""), header = T)
# 0.93
cor.test(sim.rho.50k$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.50k <- c(cor.test(sim.rho.50k$sim, unphased$sample_mean, method = "pearson")$estimate ^ 2, "mean_5")
r2.50k <- rbind(r2.50k, unphased.50k)

phased <- read.table(paste("maps/rep_1.joint.rho.50kb.bedgraph", sep = ""), header = T)
# 0.95
cor.test(sim.rho.50k$sim, phased$sample_mean, method = "spearman")$estimate
phased.50k <- c(cor.test(sim.rho.50k$sim, phased$sample_mean, method = "pearson")$estimate ^ 2, "mean_45")
r2.50k <- rbind(r2.50k, phased.50k)

mean.tmrca.vec <- numeric(length = 45)
var.tmrca.vec <- numeric(length = 45)
tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.50kb.bedgraph", sep = ""), header = T)
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
# 0.88
cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "highest.tmrca")
r2.50k <- rbind(r2.50k, high.tmrca)

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
cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "top5.tmrca")
r2.50k <- rbind(r2.50k, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.min]
}
# 0.58
cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "spearman")$estimate 
low.tmrca <- c(cor.test(sim.rho.50k$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "lowest.tmrca")
r2.50k <- rbind(r2.50k, low.tmrca)

# preparing to plot
r2.50k$r2.vec <- as.numeric(r2.50k$r2.vec)
r2.50k$type <- factor(r2.50k$type)
r2.50k$type <- factor(r2.50k$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

r2.dotplot <- ggplot(data = r2.50k, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.02, dotsize = 1, position = position_dodge(width=0.5), alpha = 0.8)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)
r2.dotplot





# 200 kb

sim.rho.200k <- read.table("../../raw_data/rig.sims.rho.2e+05.bins.txt")

nbins <- 3e+7 / 2e+5

# 45 single pairs of genomes
r2.vec <- numeric(length = 45)
rho.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".rho.200kb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.rho.200k$sim, tmp$diploid_1, method = "pearson")$estimate ^ 2
  rho.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.200k <- as.data.frame(r2.vec)
r2.200k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.rho.200kb.bedgraph", sep = ""), header = T)
# 0.93
cor.test(sim.rho.200k$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.200k <- c(cor.test(sim.rho.200k$sim, unphased$sample_mean, method = "pearson")$estimate ^ 2, "mean_5")
r2.200k <- rbind(r2.200k, unphased.200k)

phased <- read.table(paste("maps/rep_1.joint.rho.200kb.bedgraph", sep = ""), header = T)
# 0.94
cor.test(sim.rho.200k$sim, phased$sample_mean, method = "spearman")$estimate
phased.200k <- c(cor.test(sim.rho.200k$sim, unphased$sample_mean, method = "pearson")$estimate ^ 2, "mean_45")
r2.200k <- rbind(r2.200k, phased.200k)


tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.200kb.bedgraph", sep = ""), header = T)
  tmrca.df[,i] <- tmp$diploid_1
}

# builds "consensus" rec. map by cherry-picking based on local TMRCA values
pair.ids.high.tmrca <- apply(tmrca.df, 1, which.max)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.max <- pair.ids.high.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.max]
}
# 0.88
cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "highest.tmrca")
r2.200k <- rbind(r2.200k, high.tmrca)

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
cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "top5.tmrca")
r2.200k <- rbind(r2.200k, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.min]
}
# 0.67
cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.rho.200k$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "lowest.tmrca")
r2.200k <- rbind(r2.200k, low.tmrca)

# preparing to plot
r2.200k$r2.vec <- as.numeric(r2.200k$r2.vec)
r2.200k$type <- factor(r2.200k$type)
r2.200k$type <- factor(r2.200k$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

r2.dotplot <- ggplot(data = r2.200k, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.5,
                                        binwidth = 0.02, dotsize = 1, position = position_dodge(width=0.5), alpha = 0.8)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)
r2.dotplot





# 1 Mb

sim.rho.1M <- read.table("../../raw_data/rig.sims.rho.1e+06.bins.txt")

nbins <- 3e+7 / 1e+6

# 45 single pairs of genomes
r2.vec <- numeric(length = 45)
rho.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".rho.1Mb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.rho.1M$sim, tmp$diploid_1, method = "pearson")$estimate ^ 2
  rho.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.1M <- as.data.frame(r2.vec)
r2.1M$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.rho.1Mb.bedgraph", sep = ""), header = T)
# 0.89
cor.test(sim.rho.1M$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.1M <- c(cor.test(sim.rho.1M$sim, unphased$sample_mean, method = "pearson")$estimate ^ 2, "mean_5")
r2.1M <- rbind(r2.1M, unphased.1M)

phased <- read.table(paste("maps/rep_1.joint.rho.1Mb.bedgraph", sep = ""), header = T)
# 0.91
cor.test(sim.rho.1M$sim, phased$sample_mean, method = "spearman")$estimate
phased.1M <- c(cor.test(sim.rho.1M$sim, phased$sample_mean, method = "pearson")$estimate ^ 2, "mean_45")
r2.1M <- rbind(r2.1M, phased.1M)

tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.1Mb.bedgraph", sep = ""), header = T)
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
cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate
high.tmrca <- c(cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "highest.tmrca")
r2.1M <- rbind(r2.1M, high.tmrca)

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
# 0.88
cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate 
top5.tmrca <- c(cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "top5.tmrca")
r2.1M <- rbind(r2.1M, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.rec.map <- numeric(length = nrow(rho.df))
for(i in 1:nrow(rho.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.rec.map[i] <- rho.df[i, pair.min]
}
# 0.81
cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.rho.1M$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "lowest.tmrca")
r2.1M <- rbind(r2.1M, low.tmrca)

# preparing to plot
r2.1M$r2.vec <- as.numeric(r2.1M$r2.vec)
r2.1M$type <- factor(r2.1M$type)
r2.1M$type <- factor(r2.1M$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

r2.dotplot <- ggplot(data = r2.1M, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.02, dotsize = 1, position = position_dodge(width=0.5), alpha = 0.75)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)
r2.dotplot

# all together now
r2.df <- rbind(r2.50k, r2.200k, r2.1M)
r2.df$scale <- c(rep("50kb", 50), rep("200kb", 50), rep("1Mb", 50))

r2.dotplot <- ggplot(data = r2.df, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.004, dotsize = 4, position = position_dodge(width=0.3), alpha = 0.75) + ylim(0, 1)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "R2 Rec. Map")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~scale) 
ggplot2::ggsave("r2.rho.pdf", r2.dotplot, device = "pdf", width = 12, height = 9)



###################################################
#
# Theta Landscapes
#
###################################################


# 50 kb
sim.theta.50k <- read.table("../../raw_data/rig.sims.theta.50000.bins.txt")

nbins <- 3e+7 / 5e+4

# 45 single pairs of genomes
r2.vec <- numeric(length = 45)
theta.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".theta.50kb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.theta.50k$sim, tmp$diploid_1, method = "pearson")$estimate ^ 2
  theta.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.50k <- as.data.frame(r2.vec)
r2.50k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.theta.50kb.bedgraph", sep = ""), header = T)
# 0.90
cor.test(sim.theta.50k$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.50k <- c(cor.test(sim.theta.50k$sim, unphased$sample_mean, method = "pearson")$estimate ^ 2, "mean_5")
r2.50k <- rbind(r2.50k, unphased.50k)

phased <- read.table(paste("maps/rep_1.joint.theta.50kb.bedgraph", sep = ""), header = T)
# 0.90
cor.test(sim.theta.50k$sim, phased$sample_mean, method = "spearman")$estimate
phased.50k <- c(cor.test(sim.theta.50k$sim, phased$sample_mean, method = "pearson")$estimate ^ 2, "mean_45")
r2.50k <- rbind(r2.50k, phased.50k)

tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.50kb.bedgraph", sep = ""), header = T)
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
high.tmrca <- c(cor.test(sim.theta.50k$sim, cherry.pick.mut.map, method = "pearson")$estimate ^ 2, "highest.tmrca")
r2.50k <- rbind(r2.50k, high.tmrca)

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
top5.tmrca <- c(cor.test(sim.theta.50k$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "top5.tmrca")
r2.50k <- rbind(r2.50k, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.min]
}
# 0.60
cor.test(sim.theta.50k$sim, cherry.pick.mut.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.theta.50k$sim, cherry.pick.mut.map, method = "pearson")$estimate ^ 2, "lowest.tmrca")
r2.50k <- rbind(r2.50k, low.tmrca)

# preparing to plot
r2.50k$r2.vec <- as.numeric(r2.50k$r2.vec)
r2.50k$type <- factor(r2.50k$type)
r2.50k$type <- factor(r2.50k$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

r2.dotplot <- ggplot(data = r2.50k, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.02, dotsize = 1, position = position_dodge(width=0.5), alpha = 0.75)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)
r2.dotplot







# 200 kb

sim.theta.200k <- read.table("../../raw_data/rig.sims.theta.2e+05.bins.txt")

nbins <- 3e+7 / 2e+5

# 45 single pairs of genomes
r2.vec <- numeric(length = 45)
theta.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".theta.200kb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.theta.200k$sim, tmp$diploid_1, method = "pearson")$estimate ^ 2
  theta.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.200k <- as.data.frame(r2.vec)
r2.200k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.theta.200kb.bedgraph", sep = ""), header = T)
# 0.92
cor.test(sim.theta.200k$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.200k <- c(cor.test(sim.theta.200k$sim, unphased$sample_mean, method = "pearson")$estimate ^ 2, "mean_5")
r2.200k <- rbind(r2.200k, unphased.200k)

phased <- read.table(paste("maps/rep_1.joint.theta.200kb.bedgraph", sep = ""), header = T)
# 0.93
cor.test(sim.theta.200k$sim, phased$sample_mean, method = "spearman")$estimate 
phased.200k <- c(cor.test(sim.theta.200k$sim, phased$sample_mean, method = "pearson")$estimate ^ 2, "mean_45")
r2.200k <- rbind(r2.200k, phased.200k)

tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.200kb.bedgraph", sep = ""), header = T)
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
high.tmrca <- c(cor.test(sim.theta.200k$sim, cherry.pick.mut.map, method = "pearson")$estimate ^ 2, "highest.tmrca")
r2.200k <- rbind(r2.200k, high.tmrca)

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
top5.tmrca <- c(cor.test(sim.theta.200k$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "top5.tmrca")
r2.200k <- rbind(r2.200k, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.min]
}
# 0.58
cor.test(sim.theta.200k$sim, cherry.pick.mut.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.theta.200k$sim, cherry.pick.mut.map, method = "pearson")$estimate ^ 2, "lowest.tmrca")
r2.200k <- rbind(r2.200k, low.tmrca)

# preparing to plot
r2.200k$r2.vec <- as.numeric(r2.200k$r2.vec)
r2.200k$type <- factor(r2.200k$type)
r2.200k$type <- factor(r2.200k$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

r2.dotplot <- ggplot(data = r2.200k, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.02, dotsize = 1, position = position_dodge(width=0.5), alpha = 0.75)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)
r2.dotplot




# 1 Mb

sim.theta.1M <- read.table("../../raw_data/rig.sims.theta.1e+06.bins.txt")

nbins <- 3e+7 / 1e+6

# 45 single pairs of genomes
r2.vec <- numeric(length = 45)
theta.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".theta.1Mb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.theta.1M$sim, tmp$diploid_1, method = "pearson")$estimate ^ 2
  theta.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.1M <- as.data.frame(r2.vec)
r2.1M$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.theta.1Mb.bedgraph", sep = ""), header = T)
# 0.93
cor.test(sim.theta.1M$sim, unphased$sample_mean, method = "spearman")$estimate
unphased.1M <- c(cor.test(sim.theta.1M$sim, unphased$sample_mean, method = "pearson")$estimate ^ 2, "mean_5")
r2.1M <- rbind(r2.1M, unphased.1M)

phased <- read.table(paste("maps/rep_1.joint.theta.1Mb.bedgraph", sep = ""), header = T)
# 0.94
cor.test(sim.theta.1M$sim, phased$sample_mean, method = "spearman")$estimate
phased.1M <- c(cor.test(sim.theta.1M$sim, phased$sample_mean, method = "pearson")$estimate ^ 2, "mean_45")
r2.1M <- rbind(r2.1M, phased.1M)

tmrca.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.1Mb.bedgraph", sep = ""), header = T)
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
high.tmrca <- c(cor.test(sim.theta.1M$sim, cherry.pick.mut.map, method = "pearson")$estimate ^ 2, "highest.tmrca")
r2.1M <- rbind(r2.1M, high.tmrca)

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
top5.tmrca <- c(cor.test(sim.theta.1M$sim, cherry.pick.rec.map, method = "pearson")$estimate ^ 2, "top5.tmrca")
r2.1M <- rbind(r2.1M, top5.tmrca)

pair.ids.low.tmrca <- apply(tmrca.df, 1, which.min)
cherry.pick.mut.map <- numeric(length = nrow(theta.df))
for(i in 1:nrow(theta.df)) {
  pair.min <- pair.ids.low.tmrca[i]
  cherry.pick.mut.map[i] <- theta.df[i, pair.min]
}
# 0.74
cor.test(sim.theta.1M$sim, cherry.pick.mut.map, method = "spearman")$estimate
low.tmrca <- c(cor.test(sim.theta.1M$sim, cherry.pick.mut.map, method = "pearson")$estimate ^ 2, "lowest.tmrca")
r2.1M <- rbind(r2.1M, low.tmrca)

# preparing to plot
r2.1M$r2.vec <- as.numeric(r2.1M$r2.vec)
r2.1M$type <- factor(r2.1M$type)
r2.1M$type <- factor(r2.1M$type, levels = c("single_diploid", "lowest.tmrca", "highest.tmrca", "top5.tmrca", "mean_5", "mean_45"))

r2.dotplot <- ggplot(data = r2.1M, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center", stackratio = 0.6,
                                        binwidth = 0.02, dotsize = 1, position = position_dodge(width=0.5), alpha = 0.75)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)
r2.dotplot

# all together now
r2.df <- rbind(r2.50k, r2.200k, r2.1M)
r2.df$scale <- c(rep("50kb", 50), rep("200kb", 50), rep("1Mb", 50))
r2.dotplot <- ggplot(data = r2.df, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.003, dotsize = 5, position = position_dodge(width=0.3), alpha = 0.75) + ylim(0, 1)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "R2 Mut. Map")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~scale) 
ggplot2::ggsave("r2.theta.pdf", r2.dotplot, device = "pdf", width = 12, height = 9)



###################################################
#
# TMRCA Landscapes
#
###################################################

# Loads ARG:
sim.ARG <- read.tree("../../raw_data/rep_1/rep_1.ARG")

# gets nucleotide spans of each tree
arg <- readLines("../../raw_data/rep_1/rep_1.ARG")
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
r2.vec <- numeric(length = 45)
TMRCA.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.50kb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.tmrca.50k$tmrca, tmp$diploid_1, method = "pearson")$estimate ^ 2
  TMRCA.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.50k <- as.data.frame(r2.vec)
r2.50k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.TMRCA.50kb.bedgraph", sep = ""), header = T)
unphased$avg <- apply(unphased[4:ncol(unphased)], 1, mean)
cor.test(sim.tmrca.50k$tmrca, unphased$avg, method = "spearman")
unphased.50k <- c(cor.test(sim.tmrca.50k$tmrca, unphased$avg, method = "pearson")$estimate ^ 2, "mean_het_5")
r2.50k <- rbind(r2.50k, unphased.50k)

phased <- read.table(paste("maps/rep_1.joint.TMRCA.50kb.bedgraph", sep = ""), header = T)
phased$avg <- apply(phased[4:ncol(phased)], 1, mean)
cor.test(sim.tmrca.50k$tmrca, phased$avg, method = "spearman")
phased.50k <- c(cor.test(sim.tmrca.50k$tmrca, phased$avg, method = "pearson")$estimate ^ 2, "mean_het_45")
r2.50k <- rbind(r2.50k, phased.50k)


# preparing to plot
r2.50k$r2.vec <- as.numeric(r2.50k$r2.vec)
r2.50k$type <- factor(r2.50k$type)
r2.50k$type <- factor(r2.50k$type, levels = c("single_diploid", "mean_homo_5", "mean_het_5", "mean_homo_45", "mean_het_45"))

r2.dotplot <- ggplot(data = r2.50k, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.02, dotsize = 1.2, position = position_dodge(width=0.5), alpha = 0.7)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)





# 200kb
tmrca.single.nuc$bin <- ceiling(1:nrow(tmrca.single.nuc) / 200e+3)
sim.tmrca.200k <- ddply(.data = tmrca.single.nuc[-which(names(tmrca.single.nuc) == "pos")],
                        .variables = "bin", .fun = colMeans, .progress = "text")


nbins <- 3e+7 / 2e+5

# 45 single pairs of genomes
r2.vec <- numeric(length = 45)
TMRCA.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.200kb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.tmrca.200k$tmrca, tmp$diploid_1, method = "pearson")$estimate ^ 2
  TMRCA.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.200k <- as.data.frame(r2.vec)
r2.200k$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.TMRCA.200kb.bedgraph", sep = ""), header = T)
unphased$avg <- apply(unphased[4:ncol(unphased)], 1, mean)
cor.test(sim.tmrca.200k$tmrca, unphased$avg, method = "spearman")
unphased.200k <- c(cor.test(sim.tmrca.200k$tmrca, unphased$avg, method = "pearson")$estimate ^ 2, "mean_het_5")
r2.200k <- rbind(r2.200k, unphased.200k)

phased <- read.table(paste("maps/rep_1.joint.TMRCA.200kb.bedgraph", sep = ""), header = T)
phased$avg <- apply(phased[4:ncol(phased)], 1, mean)
cor.test(sim.tmrca.200k$tmrca, phased$avg, method = "spearman")
phased.200k <- c(cor.test(sim.tmrca.200k$tmrca, phased$avg, method = "pearson")$estimate ^ 2, "mean_het_45")
r2.200k <- rbind(r2.200k, phased.200k)


# preparing to plot
r2.200k$r2.vec <- as.numeric(r2.200k$r2.vec)
r2.200k$type <- factor(r2.200k$type)
r2.200k$type <- factor(r2.200k$type, levels = c("single_diploid", "mean_homo_5", "mean_het_5", "mean_homo_45", "mean_het_45"))

r2.dotplot <- ggplot(data = r2.200k, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.02, dotsize = 1.2, position = position_dodge(width=0.5), alpha = 0.7)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)







# 1Mb
tmrca.single.nuc$bin <- ceiling(1:nrow(tmrca.single.nuc) / 1e+6)
sim.tmrca.1M <- ddply(.data = tmrca.single.nuc[-which(names(tmrca.single.nuc) == "pos")],
                      .variables = "bin", .fun = colMeans, .progress = "text")

nbins <- 3e+7 / 1e+6

# 45 single pairs of genomes
r2.vec <- numeric(length = 45)
TMRCA.df <- as.data.frame(matrix(ncol = 45, nrow = nbins))
for(i in 1:45) {
  tmp <- read.table(paste("maps/rs.pair_", i, ".TMRCA.1Mb.bedgraph", sep = ""), header = T)
  r2.vec[i] <- cor.test(sim.tmrca.1M$tmrca, tmp$diploid_1, method = "pearson")$estimate ^ 2
  TMRCA.df[,i] <- tmp$diploid_1
}
summary(r2.vec)

# to organise
r2.1M <- as.data.frame(r2.vec)
r2.1M$type <- "single_diploid"

unphased <- read.table(paste("maps/rep_1.unphased.TMRCA.1Mb.bedgraph", sep = ""), header = T)
unphased$avg <- apply(unphased[4:ncol(unphased)], 1, mean)
cor.test(sim.tmrca.1M$tmrca, unphased$avg, method = "spearman")
unphased.1M <- c(cor.test(sim.tmrca.1M$tmrca, unphased$avg, method = "pearson")$estimate ^ 2, "mean_het_5")
r2.1M <- rbind(r2.1M, unphased.1M)

phased <- read.table(paste("maps/rep_1.joint.TMRCA.1Mb.bedgraph", sep = ""), header = T)
phased$avg <- apply(phased[4:ncol(phased)], 1, mean)
cor.test(sim.tmrca.1M$tmrca, phased$avg, method = "spearman")
phased.1M <- c(cor.test(sim.tmrca.1M$tmrca, phased$avg, method = "pearson")$estimate ^ 2, "mean_het_45")
r2.1M <- rbind(r2.1M, phased.1M)



# preparing to plot
r2.1M$r2.vec <- as.numeric(r2.1M$r2.vec)
r2.1M$type <- factor(r2.1M$type)
r2.1M$type <- factor(r2.1M$type, levels = c("single_diploid", "mean_homo_5", "mean_het_5", "mean_homo_45", "mean_het_45"))

r2.dotplot <- ggplot(data = r2.1M, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.02, dotsize = 1.2, position = position_dodge(width=0.5), alpha = 0.7)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "Coefficient of Determination")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 1)


# all together now
r2.df <- rbind(r2.50k, r2.200k, r2.1M)
r2.df$scale <- c(rep("50kb", 49), rep("200kb", 49), rep("1Mb", 49))
r2.dotplot <- ggplot(data = r2.df, aes(y = r2.vec, x = type, fill = type))
r2.dotplot <- r2.dotplot + geom_dotplot(binaxis = "y", binpositions = "bygroup", stackdir = "center",
                                        binwidth = 0.003, dotsize = 5, position = position_dodge(width=0.3), alpha = 0.75) + ylim(0, 1)
r2.dotplot <- r2.dotplot + labs(title = NULL, x = "Map Source", y = "R2 TMRCA Map")
r2.dotplot <- r2.dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~scale) 
ggplot2::ggsave("r2.TMRCA.pdf", r2.dotplot, device = "pdf", width = 12, height = 9)



###################################################
#
# diversity ~ rho + theta + tmrca --- 50 kb scale (UN-PHASED)
#
###################################################

no.x.axis <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
no.legend <- theme(legend.position = "none", legend.title = NULL)

# sim landscapes 50 kb
sim.rho.50k <- read.table("../../raw_data/rig.sims.rho.50000.bins.txt")
sim.theta.50k <- read.table("../../raw_data/rig.sims.theta.50000.bins.txt")

# loading computed diversity landscape
unphased.diversity.50k <- read.table("maps/rep_1.unphased.diversity.50kb.bedgraph", header = T)
unphased.diversity.50k$avg <- apply(unphased.diversity.50k[4:ncol(unphased.diversity.50k)], 1, mean)

# loading inferred landscapes
unphased.tmrca.50k <- read.table("maps/rep_1.unphased.TMRCA.50kb.bedgraph", header = T)
unphased.tmrca.50k$avg <- apply(unphased.tmrca.50k[4:ncol(unphased.tmrca.50k)], 1, mean)
unphased.rho.50k <- read.table("maps/rep_1.unphased.rho.50kb.bedgraph", header = T)
unphased.theta.50k <- read.table("maps/rep_1.unphased.theta.50kb.bedgraph", header = T)

# simulated and inferred maps are highly correlated
cor.test(x = unphased.rho.50k$sample_mean, y = sim.rho.50k$sim, method = "spearman")
cor.test(x = unphased.theta.50k$sample_mean, y = sim.theta.50k$sim, method = "spearman")

# theta and tmrca, but not rho, are positively correlated with diversity
cor.test(x = unphased.diversity.50k$avg, y = unphased.theta.50k$sample_mean, method = "spearman")
cor.test(x = unphased.diversity.50k$avg, y = unphased.tmrca.50k$avg, method = "spearman")
cor.test(x = unphased.diversity.50k$avg, y = unphased.rho.50k$sample_mean, method = "spearman")

# rho is not correlated with either theta or tmrca (good)
cor.test(x = unphased.rho.50k$sample_mean, y = unphased.theta.50k$sample_mean, method = "spearman")
cor.test(x = unphased.rho.50k$sample_mean, y = unphased.tmrca.50k$avg, method = "spearman")
# theta is still weakly correlated with tmrca 
cor.test(x = unphased.theta.50k$sample_mean, y = unphased.tmrca.50k$avg, method = "spearman")


inf.lands <- as.data.frame(cbind(unphased.diversity.50k$avg, unphased.theta.50k$sample_mean, unphased.rho.50k$sample_mean, unphased.tmrca.50k$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)



# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 50, y = value, colour = "#F8766D")) 
diversity.map <- diversity.map + geom_line(data = molten.diversity)
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D")
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi)) + no.legend
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, 2 * inf.lands$rho, sim.rho.50k$sim)) # adds simulated rec. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 50, y = value, colour = variable)) 
rho.map <- rho.map + geom_line(data = molten.rho) + scale_fill_manual(values = c("#fc8d59", "#99d594"))
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.009, label = "R_2 = 67.8%")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.50k$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 50, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.01, label = "R_2 = 76.4%")

tmrca.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$tmrca, sim.tmrca.50k$tmrca))
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
lands.plot <-  plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::ggsave("landscapes.50kb.unphased.pdf", lands.plot, width = 9, height = 15)







# Linear model WITHOUT interactions

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

anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta       1 0.66540 0.66540 15039.111 0.000000 0.75362
# rho         1 0.00044 0.00044    10.053 0.001599 0.00050
# tmrca       1 0.19073 0.19073  4310.741 0.000000 0.21601
# Residuals 596 0.02637 0.00004                    0.02987

relimp.diversity <- calc.relimp(m.diversity.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 0.001474016 
# Proportion of variance explained by model: 97.01%
# Relative importance metrics: 
# theta 0.64273197
# rho   0.00522212
# tmrca 0.35204591






# Linear model with interactions

m.diversity.int <- lm(diversity ~ (theta + rho + tmrca) ^ 2, data = inf.lands)
plot(m.diversity.int, which = 2)

bc.diversity <- boxcox(m.diversity.int, lambda = seq(0.5, 1.5, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.int.bc <- update(m.diversity.int, (diversity^l -1)/l~.)
plot(m.diversity.int.bc, which = 2)

shapiro.test(resid(m.diversity.int.bc)) #***
hist(resid(m.diversity.int.bc), nclass = 30) # not too bad?
hmctest(m.diversity.int.bc, nsim = 3000) # **
dwtest(m.diversity.int.bc) # ***
Box.test(resid(m.diversity.int.bc)[order(predict(m.diversity.int.bc))], type = "Ljung-Box") # NS

summary(m.diversity.int.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.990e+00  5.197e-03 -382.881  < 2e-16 ***
# theta        2.993e-01  4.375e-02     6.841 1.96e-11 ***
# rho         -3.354e-01  2.152e-01    -1.559   0.1196    
# tmrca        1.675e-03  3.286e-04     5.096 4.67e-07 ***
# theta:rho   -3.847e+01  1.543e+01    -2.493   0.0129 *  
# theta:tmrca  2.026e+00  5.267e-02    38.454  < 2e-16 ***
# rho:tmrca    1.134e+00  2.854e-01     3.972 8.00e-05 ***
# Multiple R-squared:  0.9866,	Adjusted R-squared:  0.9864 

anova.diversity <- anova(m.diversity.int.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta         1 0.0140245 0.0140245 32298.630 0.0000e+00 0.73242
# rho           1 0.0000462 0.0000462   106.334 0.0000e+00 0.00241
# tmrca         1 0.0041493 0.0041493  9555.922 0.0000e+00 0.21669
# theta:rho     1 0.0000182 0.0000182    41.801 0.0000e+00 0.00095
# theta:tmrca   1 0.0006457 0.0006457  1487.042 0.0000e+00 0.03372
# rho:tmrca     1 0.0000069 0.0000069    15.777 8.0028e-05 0.00036
# Residuals   593 0.0002575 0.0000004                      0.01345

relimp.diversity <- calc.relimp(m.diversity.int.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 3.196695e-05 
# Proportion of variance explained by model: 98.66%
# Relative importance metrics: 
# theta       0.6149844232
# rho         0.0029941098
# tmrca       0.3456617073
# theta:rho   0.0014856720
# theta:tmrca 0.0344348214
# rho:tmrca   0.0004392663


# Comparing the model with interactions and the model without
anova(m.diversity.bc, m.diversity.int.bc) #***


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
sim.rho.200k <- read.table("../../raw_data/rig.sims.rho.2e+05.bins.txt")
sim.theta.200k <- read.table("../../raw_data/rig.sims.theta.2e+05.bins.txt")

# loading computed diversity landscape
unphased.diversity.200k <- read.table("maps/rep_1.unphased.diversity.200kb.bedgraph", header = T)
unphased.diversity.200k$avg <- apply(unphased.diversity.200k[4:ncol(unphased.diversity.200k)], 1, mean)

# loading inferred landscapes
unphased.tmrca.200k <- read.table("maps/rep_1.unphased.TMRCA.200kb.bedgraph", header = T)
unphased.tmrca.200k$avg <- apply(unphased.tmrca.200k[4:ncol(unphased.tmrca.200k)], 1, mean)
unphased.rho.200k <- read.table("maps/rep_1.unphased.rho.200kb.bedgraph", header = T)
unphased.theta.200k <- read.table("maps/rep_1.unphased.theta.200kb.bedgraph", header = T)

# simulated and inferred maps are highly correlated
cor.test(x = unphased.rho.200k$sample_mean, y = sim.rho.200k$sim, method = "spearman")
cor.test(x = unphased.theta.200k$sample_mean, y = sim.theta.200k$sim, method = "spearman")

# theta and tmrca, but not rho, are positively correlated with diversity
cor.test(x = unphased.diversity.200k$avg, y = unphased.theta.200k$sample_mean, method = "spearman")
cor.test(x = unphased.diversity.200k$avg, y = unphased.tmrca.200k$avg, method = "spearman")
cor.test(x = unphased.diversity.200k$avg, y = unphased.rho.200k$sample_mean, method = "spearman")

# rho is not correlated with either theta or tmrca (good)
cor.test(x = unphased.rho.200k$sample_mean, y = unphased.theta.200k$sample_mean, method = "spearman")
cor.test(x = unphased.rho.200k$sample_mean, y = unphased.tmrca.200k$avg, method = "spearman")
# theta is weakly correlated with tmrca (not so good)
cor.test(x = unphased.theta.200k$sample_mean, y = unphased.tmrca.200k$avg, method = "spearman")

inf.lands <- as.data.frame(cbind(unphased.diversity.200k$avg, unphased.theta.200k$sample_mean, unphased.rho.200k$sample_mean, unphased.tmrca.200k$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

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
rho.map <- rho.map + annotate("text", x = 15000, y = 0.0075, label = "R_2 = 74.4%")

theta.plot <-  as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.200k$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 200, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta)
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.01, label = "R_2 = 79.4%")

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
lands.plot <-  plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::ggsave("landscapes.200kb.unphased.pdf", lands.plot, width = 9, height = 15)







# Linear model WITHOUT interactions

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands)
plot(m.diversity, which = 2)

shapiro.test(resid(m.diversity)) #**
hist(resid(m.diversity), nclass = 30) # not too bad?
hmctest(m.diversity, nsim = 3000) # NS
dwtest(m.diversity) # NS
Box.test(resid(m.diversity)[order(predict(m.diversity))], type = "Ljung-Box") #

summary(m.diversity)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.0039790  0.0001911 -20.824   <2e-16 ***
# theta        0.7315314  0.0180289  40.575   <2e-16 ***
# rho         -0.0564721  0.0622994  -0.906    0.366    
# tmrca        0.0055288  0.0002552  21.662   <2e-16 ***
# Multiple R-squared:  0.9514,	Adjusted R-squared:  0.9504

anova.diversity <- anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta       1 0.00048991 0.00048991 2377.45 0.00000000 0.79112
# rho         1 0.00000257 0.00000257   12.48 0.00055104 0.00415
# tmrca       1 0.00009669 0.00009669  469.23 0.00000000 0.15614
# Residuals 146 0.00003009 0.00000021                    0.04858

relimp.diversity <- calc.relimp(m.diversity, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 0.007611425 
# Proportion of variance explained by model: 95.14%
# Relative importance metrics: 
# theta 0.705773765
# rho   0.004109551
# tmrca 0.290116684






# Linear model with interactions

m.diversity.int <- lm(diversity ~ (theta + rho + tmrca) ^ 2, data = inf.lands)
plot(m.diversity.int, which = 2)

bc.diversity <- boxcox(m.diversity.int, lambda = seq(-0.5, 0.5, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.int.bc <- update(m.diversity.int, (diversity^l -1)/l~.)
plot(m.diversity.int.bc, which = 2)

shapiro.test(resid(m.diversity.int.bc)) #*
hist(resid(m.diversity.int.bc), nclass = 30) # not too bad?
hmctest(m.diversity.int.bc, nsim = 3000) # NS
dwtest(m.diversity.int.bc) # **
Box.test(resid(m.diversity.int.bc)[order(predict(m.diversity.int.bc))], type = "Ljung-Box") # *

summary(m.diversity.int.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.990e+00  5.197e-03 -382.881  < 2e-16 ***
# theta        9.301e+00  9.277e-01   10.026  < 2e-16 ***
# rho         -6.628e+00  5.069e+00   -1.308  0.19309    
# tmrca        5.428e-02  7.092e-03    7.654 2.64e-12 ***
# theta:rho   -7.722e+02  2.965e+02   -2.604  0.01018 *  
# theta:tmrca  4.862e+00  1.152e+00    4.221 4.30e-05 ***
# rho:tmrca    1.847e+01  6.763e+00    2.730  0.00712 ** 
# Multiple R-squared:  0.9798,	Adjusted R-squared:  0.979

anova.diversity <- anova(m.diversity.int.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta         1 0.135587 0.135587 5800.4180 0.0000000 0.81768
# rho           1 0.000161 0.000161    6.8964 0.0095773 0.00097
# tmrca         1 0.025769 0.025769 1102.4083 0.0000000 0.15541
# theta:rho     1 0.000315 0.000315   13.4850 0.0003389 0.00190
# theta:tmrca   1 0.000469 0.000469   20.0728 0.0000152 0.00283
# rho:tmrca     1 0.000174 0.000174    7.4543 0.0071248 0.00105
# Residuals   143 0.003343 0.000023                     0.02016

relimp.diversity <- calc.relimp(m.diversity.int.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 0.007611425 
# Proportion of variance explained by model: 97.98%
# Relative importance metrics: 
# theta       0.706570681
# rho         0.006820388
# tmrca       0.279731933
# theta:rho   0.002601175
# theta:tmrca 0.003003132
# rho:tmrca   0.001272692


# Comparing the model with interactions and the model without
anova(m.diversity, m.diversity.int) #***

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
sim.rho.1M <- read.table("../../raw_data/rig.sims.rho.1e+06.bins.txt")
sim.theta.1M <- read.table("../../raw_data/rig.sims.theta.1e+06.bins.txt")

# loading computed diversity landscape
unphased.diversity.1M <- read.table("maps/rep_1.unphased.diversity.1Mb.bedgraph", header = T)
unphased.diversity.1M$avg <- apply(unphased.diversity.1M[4:ncol(unphased.diversity.1M)], 1, mean)

# loading inferred landscapes
unphased.tmrca.1M <- read.table("maps/rep_1.unphased.TMRCA.1Mb.bedgraph", header = T)
unphased.tmrca.1M$avg <- apply(unphased.tmrca.1M[4:ncol(unphased.tmrca.1M)], 1, mean)
unphased.rho.1M <- read.table("maps/rep_1.unphased.rho.1Mb.bedgraph", header = T)
unphased.theta.1M <- read.table("maps/rep_1.unphased.theta.1Mb.bedgraph", header = T)

# simulated and inferred maps are highly correlated
cor.test(x = unphased.rho.1M$sample_mean, y = sim.rho.1M$sim, method = "spearman")
cor.test(x = unphased.theta.1M$sample_mean, y = sim.theta.1M$sim, method = "spearman")

# theta and tmrca, but not rho, are positively correlated with diversity
cor.test(x = unphased.diversity.1M$avg, y = unphased.theta.1M$sample_mean, method = "spearman")
cor.test(x = unphased.diversity.1M$avg, y = unphased.tmrca.1M$avg, method = "spearman")
cor.test(x = unphased.diversity.1M$avg, y = unphased.rho.1M$sample_mean, method = "spearman")

# rho is not correlated with either theta or tmrca (good)
cor.test(x = unphased.rho.1M$sample_mean, y = unphased.theta.1M$sample_mean, method = "spearman")
cor.test(x = unphased.rho.1M$sample_mean, y = unphased.tmrca.1M$avg, method = "spearman")
# at the 1Mb scale, theta is NOT correlated with tmrca (good)
cor.test(x = unphased.theta.1M$sample_mean, y = unphased.tmrca.1M$avg, method = "spearman")

inf.lands <- as.data.frame(cbind(unphased.diversity.1M$avg, unphased.theta.1M$sample_mean, unphased.rho.1M$sample_mean, unphased.tmrca.1M$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

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
rho.map <- rho.map + annotate("text", x = 15000, y = 0.0025, label = "R_2 = 78.6%")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.1M$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 1000, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta)
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.0085, label = "R_2 = 84.8%")

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
lands.plot <-  plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
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
r2.tab <- as.data.frame(matrix(ncol = 5, nrow = 6))
colnames(r2.tab) <- c("Total", "Theta", "Rho", "TMRCA", "Bin_size(kb)")




# sim landscapes 50 kb
sim.rho.50k <- read.table("../../raw_data/rig.sims.rho.50000.bins.txt")
sim.theta.50k <- read.table("../../raw_data/rig.sims.theta.50000.bins.txt")

# loading computed diversity landscape
joint.diversity.50k <- read.table("maps/rep_1.joint.diversity.50kb.bedgraph", header = T)
joint.diversity.50k$avg <- apply(joint.diversity.50k[4:ncol(joint.diversity.50k)], 1, mean)

# loading inferred landscapes
joint.tmrca.50k <- read.table("maps/rep_1.joint.TMRCA.50kb.bedgraph", header = T)
joint.tmrca.50k$avg <- apply(joint.tmrca.50k[4:ncol(joint.tmrca.50k)], 1, mean)
joint.rho.50k <- read.table("maps/rep_1.joint.rho.50kb.bedgraph", header = T)
joint.theta.50k <- read.table("maps/rep_1.joint.theta.50kb.bedgraph", header = T)

# theta and tmrca, but not rho, are positively correlated with diversity
cor.test(x = joint.diversity.50k$avg, y = joint.theta.50k$sample_mean, method = "spearman")
cor.test(x = joint.diversity.50k$avg, y = joint.tmrca.50k$avg, method = "spearman")
cor.test(x = joint.diversity.50k$avg, y = joint.rho.50k$sample_mean, method = "spearman")

# rho is not correlated with either theta or tmrca (good)
cor.test(x = joint.rho.50k$sample_mean, y = joint.theta.50k$sample_mean, method = "spearman")
cor.test(x = joint.rho.50k$sample_mean, y = joint.tmrca.50k$avg, method = "spearman")
# theta is correlated with tmrca 
cor.test(x = joint.theta.50k$sample_mean, y = joint.tmrca.50k$avg, method = "spearman")

inf.lands <- as.data.frame(cbind(joint.diversity.50k$avg, joint.theta.50k$sample_mean, joint.rho.50k$sample_mean, joint.tmrca.50k$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)



# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 50, y = value, colour = "#F8766D")) 
diversity.map <- diversity.map + geom_line(data = molten.diversity)
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, 2 * inf.lands$rho, sim.rho.50k$sim)) # adds simulated rec. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 50, y = value, colour = variable)) 
rho.map <- rho.map + geom_line(data = molten.rho) + scale_fill_manual(values = c("#fc8d59", "#99d594"))
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.009, label = "R_2 = 69.3%")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.50k$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 50, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta) + scale_fill_manual(values = c("#fee08b", "#3288bd"))
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.0095, label = "R_2 = 79.9%")

tmrca.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$tmrca, sim.tmrca.50k$tmrca))
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
lands.plot <-  plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::ggsave("landscapes.50kb.phased.pdf", lands.plot, width = 9, height = 15)







# Linear model WITHOUT interactions

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

r2.tab[1,] <- c(97.62, 76.12, 0, 21.46, 50)

anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta       1 0.48113 0.48113 19094.739 0.00000000 0.76119
# rho         1 0.00030 0.00030    11.824 0.00062561 0.00047
# tmrca       1 0.13563 0.13563  5382.881 0.00000000 0.21458

relimp.diversity <- calc.relimp(m.diversity.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 0.001055215 
# Proportion of variance explained by model: 97.62%
# Relative importance metrics: 
# theta       1 0.48113 0.48113 19094.739 0.00000000 0.76119
# rho         1 0.00030 0.00030    11.824 0.00062561 0.00047
# tmrca       1 0.13563 0.13563  5382.881 0.00000000 0.21458
# Residuals 596 0.01502 0.00003                      0.02376


# Linear model with interactions

m.diversity.int <- lm(diversity ~ (theta + rho + tmrca) ^ 2, data = inf.lands)
plot(m.diversity.int, which = 2)

bc.diversity <- boxcox(m.diversity.int, lambda = seq(0.5, 1.5, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.int.bc <- update(m.diversity.int, (diversity^l -1)/l~.)
plot(m.diversity.int.bc, which = 2)

shapiro.test(resid(m.diversity.int.bc)) #***
hist(resid(m.diversity.int.bc), nclass = 30) # not too bad?
hmctest(m.diversity.int.bc, nsim = 3000) # *
dwtest(m.diversity.int.bc) # ***
Box.test(resid(m.diversity.int.bc)[order(predict(m.diversity.int.bc))], type = "Ljung-Box") # NS

summary(m.diversity.int.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.232e+00  2.755e-04 -4470.326  < 2e-16 ***
# theta        5.114e-01  7.194e-02     7.108 3.39e-12 ***
# rho          2.339e-01  3.356e-01     0.697 0.486092    
# tmrca        2.453e-03  2.982e-04     8.225 1.24e-15 ***
# theta:rho   -1.030e+02  3.008e+01    -3.425 0.000656 ***
# theta:tmrca  2.482e+00  6.979e-02    35.567  < 2e-16 ***
# rho:tmrca    5.363e-01  3.440e-01     1.559 0.119584    
# Multiple R-squared: 0.9886,	Adjusted R-squared:  0.9885 

anova.diversity <- anova(m.diversity.int.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta         1 0.0182731 0.0182731 39817.9871 0.00000 0.76484
# rho           1 0.0000001 0.0000001     0.2908 0.58989 0.00001
# tmrca         1 0.0047347 0.0047347 10317.0933 0.00000 0.19817
# theta:rho     1 0.0000173 0.0000173    37.6436 0.00000 0.00072
# theta:tmrca   1 0.0005930 0.0005930  1292.0941 0.00000 0.02482
# rho:tmrca     1 0.0000011 0.0000011     2.4298 0.11958 0.00005
# Residuals   593 0.0002721 0.0000005                    0.01139

relimp.diversity <- calc.relimp(m.diversity.int.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance:  3.988547e-05 
# Proportion of variance explained by model: 98.86%
# Relative importance metrics: 
# theta       0.6331686096
# rho         0.0030113073
# tmrca       0.3367418685
# theta:rho   0.0017053132
# theta:tmrca 0.0250433273
# rho:tmrca   0.0003295742


# Comparing the model with interactions and the model without
anova(m.diversity.bc, m.diversity.int.bc) #***


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
sim.rho.200k <- read.table("../../raw_data/rig.sims.rho.2e+05.bins.txt")
sim.theta.200k <- read.table("../../raw_data/rig.sims.theta.2e+05.bins.txt")

# loading computed diversity landscape
joint.diversity.200k <- read.table("maps/rep_1.joint.diversity.200kb.bedgraph", header = T)
joint.diversity.200k$avg <- apply(joint.diversity.200k[4:ncol(joint.diversity.200k)], 1, mean)

# loading inferred landscapes
joint.tmrca.200k <- read.table("maps/rep_1.joint.TMRCA.200kb.bedgraph", header = T)
joint.tmrca.200k$avg <- apply(joint.tmrca.200k[4:ncol(joint.tmrca.200k)], 1, mean)
joint.rho.200k <- read.table("maps/rep_1.joint.rho.200kb.bedgraph", header = T)
joint.theta.200k <- read.table("maps/rep_1.joint.theta.200kb.bedgraph", header = T)

# theta and tmrca, but not rho, are positively correlated with diversity
cor.test(x = joint.diversity.200k$avg, y = joint.theta.200k$sample_mean, method = "spearman")
cor.test(x = joint.diversity.200k$avg, y = joint.tmrca.200k$avg, method = "spearman")
cor.test(x = joint.diversity.200k$avg, y = joint.rho.200k$sample_mean, method = "spearman")

# rho is not correlated with either theta or tmrca (good)
cor.test(x = joint.rho.200k$sample_mean, y = joint.theta.200k$sample_mean, method = "spearman")
cor.test(x = joint.rho.200k$sample_mean, y = joint.tmrca.200k$avg, method = "spearman")
# theta is weakly correlated with tmrca (not so good)
cor.test(x = joint.theta.200k$sample_mean, y = joint.tmrca.200k$avg, method = "spearman")

inf.lands <- as.data.frame(cbind(joint.diversity.200k$avg, joint.theta.200k$sample_mean, joint.rho.200k$sample_mean, joint.tmrca.200k$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 200, y = value, colour = "#F8766D")) 
diversity.map <- diversity.map + geom_line(data = molten.diversity)
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, 2 * inf.lands$rho, sim.rho.200k$sim)) # adds simulated rec. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 200, y = value, colour = variable)) 
rho.map <- rho.map + geom_line(data = molten.rho)
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.007, label = "R_2 = 74.4%")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.200k$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 200, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta)
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.0085, label = "R_2 = 83.9%")

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
lands.plot <-  plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::ggsave("landscapes.200kb.phased.pdf", lands.plot, width = 9, height = 15)







# Linear model WITHOUT interactions

m.diversity <- lm(diversity ~ theta + rho + tmrca, data = inf.lands)
plot(m.diversity, which = 2)

bc.diversity <- boxcox(m.diversity, lambda = seq(0., 1., len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.bc <- update(m.diversity, (diversity^l -1)/l~.)
plot(m.diversity.int.bc, which = 2)

shapiro.test(resid(m.diversity.bc)) #*
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

r2.tab[2,] <- c(97.97, 83.94, 0, 14.1, 200)

anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df   Sum Sq  Mean Sq   F value  Pr(>F)  VarExp
# theta       1 0.059395 0.059395 6153.9711 0.00000 0.83942
# rho         1 0.000007 0.000007    0.7615 0.38428 0.00010
# tmrca       1 0.009946 0.009946 1030.5165 0.00000 0.14056
# Residuals 146 0.001409 0.000010                   0.01991

relimp.diversity <- calc.relimp(m.diversity.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 0.0004748805 
# Proportion of variance explained by model: 98.01%
# Relative importance metrics: 
# theta 0.7243439
# rho   0.0042625
# tmrca 0.2713936






# Linear model with interactions

m.diversity.int <- lm(diversity ~ (theta + rho + tmrca) ^ 2, data = inf.lands)
plot(m.diversity.int, which = 2)

bc.diversity <- boxcox(m.diversity.int, lambda = seq(0.5, 1.5, len = 500))
l <- bc.diversity$x[which.max(bc.diversity$y)]
m.diversity.int.bc <- update(m.diversity.int, (diversity^l -1)/l~.)
plot(m.diversity.int.bc, which = 2)

shapiro.test(resid(m.diversity.int.bc)) #**
hist(resid(m.diversity.int.bc), nclass = 30) # not too bad?
hmctest(m.diversity.int.bc, nsim = 3000) # NS
dwtest(m.diversity.int.bc) # *
Box.test(resid(m.diversity.int.bc)[order(predict(m.diversity.int.bc))], type = "Ljung-Box") # NS

summary(m.diversity.int.bc)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.990e+00  5.197e-03 -382.881  < 2e-16 ***
# theta        5.833e-01  1.511e-01     3.861 0.000171 ***
# rho          5.142e-01  6.962e-01     0.739 0.461380    
# tmrca        1.884e-03  6.348e-04     2.968 0.003517 ** 
# theta:rho   -1.533e+02  5.750e+01    -2.666 0.008564 ** 
# theta:tmrca  1.996e+00  1.495e-01    13.353  < 2e-16 ***
# rho:tmrca    2.655e-01  7.160e-01     0.371 0.711340 
# Multiple R-squared: 0.9887,	Adjusted R-squared:  0.9882 

anova.diversity <- anova(m.diversity.int.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta         1 0.00251842 0.00251842 10569.7628 0.00000 0.83639
# rho           1 0.00000015 0.00000015     0.6487 0.42190 0.00005
# tmrca         1 0.00040896 0.00040896  1716.3962 0.00000 0.13582
# theta:rho     1 0.00000572 0.00000572    24.0273 0.00000 0.00190
# theta:tmrca   1 0.00004370 0.00004370   183.3909 0.00000 0.01451
# rho:tmrca     1 0.00000003 0.00000003     0.1375 0.71134 0.00001
# Residuals   143 0.00003407 0.00000024                    0.01132

relimp.diversity <- calc.relimp(m.diversity.int.bc, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 2.020849e-05 
# Proportion of variance explained by model: 98.87%
# Relative importance metrics: 
# theta       0.7175773764
# rho         0.0029034912
# tmrca       0.2605879858
# theta:rho   0.0029241610
# theta:tmrca 0.0151223840
# rho:tmrca   0.0008846016

# Comparing the model with interactions and the model without
anova(m.diversity.bc, m.diversity.int.bc) #***


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
sim.rho.1M <- read.table("../../raw_data/rig.sims.rho.1e+06.bins.txt")
sim.theta.1M <- read.table("../../raw_data/rig.sims.theta.1e+06.bins.txt")

# loading computed diversity landscape
joint.diversity.1M <- read.table("maps/rep_1.joint.diversity.1Mb.bedgraph", header = T)
joint.diversity.1M$avg <- apply(joint.diversity.1M[4:ncol(joint.diversity.1M)], 1, mean)

# loading inferred landscapes
joint.tmrca.1M <- read.table("maps/rep_1.joint.TMRCA.1Mb.bedgraph", header = T)
joint.tmrca.1M$avg <- apply(joint.tmrca.1M[4:ncol(joint.tmrca.1M)], 1, mean)
joint.rho.1M <- read.table("maps/rep_1.joint.rho.1Mb.bedgraph", header = T)
joint.theta.1M <- read.table("maps/rep_1.joint.theta.1Mb.bedgraph", header = T)

# theta and tmrca, but not rho, are positively correlated with diversity
cor.test(x = joint.diversity.1M$avg, y = joint.theta.1M$sample_mean, method = "spearman")
cor.test(x = joint.diversity.1M$avg, y = joint.tmrca.1M$avg, method = "spearman")
cor.test(x = joint.diversity.1M$avg, y = joint.rho.1M$sample_mean, method = "spearman")

# rho is not correlated with either theta or tmrca (good)
cor.test(x = joint.rho.1M$sample_mean, y = joint.theta.1M$sample_mean, method = "spearman")
cor.test(x = joint.rho.1M$sample_mean, y = joint.tmrca.1M$avg, method = "spearman")
# at the 1Mb scale, theta is NOT correlated with tmrca (good)
cor.test(x = joint.theta.1M$sample_mean, y = joint.tmrca.1M$avg, method = "spearman")

inf.lands <- as.data.frame(cbind(joint.diversity.1M$avg, joint.theta.1M$sample_mean, joint.rho.1M$sample_mean, joint.tmrca.1M$avg))
names(inf.lands) <- c("diversity", "theta", "rho", "tmrca")
inf.lands$bin <- 1:nrow(inf.lands)

# plots
scale.3d <- function(x) sprintf("%.3f", x) # digits shown in y axis

molten.diversity <- melt(inf.lands[c(1,5)], id.vars = "bin")
diversity.map <- ggplot(data = molten.diversity, aes(x = bin * 1000, y = value)) 
diversity.map <- diversity.map + geom_line(data = molten.diversity, colour = "#F8766D")
diversity.map <- diversity.map + geom_smooth(method = "loess", se = F, colour = "#F8766D") + no.legend
diversity.map <- diversity.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
diversity.map <- diversity.map + labs(title = NULL, x = NULL, y = expression(pi))
diversity.map <- diversity.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

rho.plot <- as.data.frame(cbind(inf.lands$bin, 2 * inf.lands$rho, sim.rho.1M$sim)) # adds simulated mut. map
names(rho.plot) <- c("bin", "inf", "sim")
molten.rho <- melt(rho.plot, id.vars = "bin")
rho.map <- ggplot(data = molten.rho, aes(x = bin * 1000, y = value, colour = variable)) 
rho.map <- rho.map + geom_line(data = molten.rho)
rho.map <- rho.map + geom_smooth(method = "loess", se = F) + no.legend
rho.map <- rho.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
rho.map <- rho.map + labs(title = NULL, x = NULL, y = expression(rho))
rho.map <- rho.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
rho.map <- rho.map + annotate("text", x = 15000, y = 0.0025, label = "R_2 = 77.8%")

theta.plot <- as.data.frame(cbind(inf.lands$bin, inf.lands$theta, sim.theta.1M$sim)) # adds simulated mut. map
names(theta.plot) <- c("bin", "inf", "sim")
molten.theta <- melt(theta.plot, id.vars = "bin")
theta.map <- ggplot(data = molten.theta, aes(x = bin * 1000, y = value, colour = variable)) 
theta.map <- theta.map + geom_line(data = molten.theta)
theta.map <- theta.map + geom_smooth(method = "loess", se = F) + no.legend
theta.map <- theta.map + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks(), labels = scale.3d)
theta.map <- theta.map + labs(title = NULL, x = NULL, y = expression(theta))
theta.map <- theta.map + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
theta.map <- theta.map + annotate("text", x = 15000, y = 0.006, label = "R_2 = 89.7%")

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
lands.plot <-  plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.1), scale = c(1, 0.2))
cowplot::ggsave("landscapes.1Mb.phased.pdf", lands.plot, width = 9, height = 15)







# Linear model WITHOUT interactions

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

r2.tab[3,] <- c(96.54, 89.1, 0.6, 7.16, 1000)

anova.diversity <- anova(m.diversity)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df     Sum Sq    Mean Sq  F value    Pr(>F)  VarExp
# theta      1 3.3265e-05 3.3265e-05 746.397 0.000000 0.89076
# rho        1 2.4700e-07 2.4700e-07   5.544 0.026375 0.00662
# tmrca      1 2.6740e-06 2.6740e-06  59.994 0.000000 0.07160
# Residuals 26 1.1590e-06 4.5000e-08                  0.03103

relimp.diversity <- calc.relimp(m.diversity, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 1.287749e-06 
# Proportion of variance explained by model: 96.9%
# Relative importance metrics: 
# theta 0.84478752
# rho   0.01109592
# tmrca 0.14411656






# Linear model with interactions

m.diversity.int <- lm(diversity ~ (theta + rho + tmrca) ^ 2, data = inf.lands)
plot(m.diversity.int, which = 2)

shapiro.test(resid(m.diversity.int)) #*
hist(resid(m.diversity.int.bc), nclass = 30)
hmctest(m.diversity.int.bc, nsim = 3000) # NS
dwtest(m.diversity.int.bc) # NS
Box.test(resid(m.diversity.int.bc)[order(predict(m.diversity.int.bc))], type = "Ljung-Box") # NS

summary(m.diversity.int)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.002698   0.001658   1.628 0.117185    
# theta        -1.111257   0.552304  -2.012 0.056072 .  
# rho           0.990009   1.748992   0.566 0.576844    
# tmrca        -0.003096   0.001737  -1.782 0.087886 .  
# theta:rho     6.318193 117.986814   0.054 0.957756    
# theta:tmrca   2.273322   0.578728   3.928 0.000672 ***
# rho:tmrca    -1.071241   1.755685  -0.610 0.547741    
# Multiple R-squared: 0.9834,	Adjusted R-squared:  0.9791 

anova.diversity <- anova(m.diversity.int)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)

anova.diversity
# Analysis of Variance Table
# Df  Sum Sq Mean Sq F value     Pr(>F)  VarExp
# theta        1 3.3265e-05 3.3265e-05 1233.1206 0.00000 0.89076
# rho          1 2.4700e-07 2.4700e-07    9.1592 0.00601 0.00662
# tmrca        1 2.6740e-06 2.6740e-06   99.1153 0.00000 0.07160
# theta:rho    1 8.4000e-08 8.4000e-08    3.1149 0.09086 0.00225
# theta:tmrca  1 4.4400e-07 4.4400e-07   16.4673 0.00049 0.01190
# rho:tmrca    1 1.0000e-08 1.0000e-08    0.3723 0.54774 0.00027
# Residuals   23 6.2000e-07 2.7000e-08                   0.01661

relimp.diversity <- calc.relimp(m.diversity.int, rela = TRUE, type = "lmg")

relimp.diversity
# Total response variance: 1.287749e-06 
# Proportion of variance explained by model: 98.34%
# Relative importance metrics: 
# theta       0.830691671
# rho         0.010785464
# tmrca       0.140702021
# theta:rho   0.002191028
# theta:tmrca 0.012727245
# rho:tmrca   0.0029025726


# Comparing the model with interactions and the model without
anova(m.diversity, m.diversity.int) #**





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
anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)
anova.diversity
# theta       1 0.110534 0.110534 7410.5194 0.000000 0.77100
# rho         1 0.000042 0.000042    2.8372 0.092631 0.00030
# tmrca       1 0.023898 0.023898 1602.1907 0.000000 0.16669
# Residuals 596 0.008890 0.000015                    0.06201

r2.tab[4,] <- c(93.77, 77.1, 0, 16.70, 50)



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
anova.diversity <- anova(m.diversity.bc)
apiss <- anova.diversity$"Sum Sq"
anova.diversity$VarExp <- apiss / sum(apiss)
anova.diversity
# theta       1 0.0065549 0.0065549 3129.142 0.000000 0.85259
# rho         1 0.0000123 0.0000123    5.873 0.016598 0.00160
# tmrca       1 0.0008152 0.0008152  389.166 0.000000 0.10603
# Residuals 146 0.0003058 0.0000021                   0.03978

r2.tab[5,] <- c(95.94, 85.3, 0, 10.60, 200)





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

r2.tab[6,] <- c(96.6, 91.93, 0, 4.90, 1000)


########################################
#
# R2 Plot
#
########################################

r2.tab.2 <- as.data.frame(cbind(apply(r2.tab, 2, as.numeric)))
names(r2.tab.2)[5] <- "bin.size"
r2.tab.2$type <- c(rep("inf", 3), rep("sim", 3))

molten.r2 <- melt(r2.tab.2, id.vars = c("bin.size", "type"))
r2.plot <- ggplot(data = molten.r2, aes(x = bin.size, y = value, colour = variable, fill = type))
r2.plot <- r2.plot + geom_line(data = molten.r2)
r2.plot <- r2.plot + geom_point(aes(shape = type, colour = variable), size = 5)
r2.plot <- r2.plot + scale_x_continuous(breaks = c(50, 200, 1000))
r2.plot <- r2.plot + scale_y_continuous(breaks = pretty_breaks())
r2.plot <- r2.plot + labs(title = NULL, x = "Bin Size (kb)", y = "Var. Explained")
r2.plot <- r2.plot + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

ggsave("lm.r2.sim.pdf", r2.plot, width = 10, height = 10)

