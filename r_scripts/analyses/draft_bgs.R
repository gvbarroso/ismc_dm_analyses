
#Regular Region

# true landscapes, 1r

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

# true landscapes, 10r

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

# true landscapes, 1r

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

# true landscapes, 10r

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
r2.true.avg$r <- as.factor(rep(c(rep("1r", 3), rep("10r", 3)), 2))
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

