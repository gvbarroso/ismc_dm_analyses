# Created: 17/06/2019
# Last modified: 17/06/2019
# Author: Gustavo Barroso 
# This script bins simulated theta and rho landscapes 

library(plyr)

# function to read simulated landscape in coordinate/fragments and return a rho value per position
expand_sim_landscape <- function(sim_landscape, sequence_length) {
  
  sim <- vector(length = sequence_length)
  # makes continuous landscape
  num_of_regions <- nrow(sim_landscape)
  for(i in 1:num_of_regions){
    if(i == 1) {
      region_end <- sim_landscape$V2[2]
      sim[1:region_end] <- sim_landscape$V1[1]
    } 
    else if(i == num_of_regions) {
      region_start <- sim_landscape$V2[i]
      region_end <- sequence_length
      sim[region_start:region_end] <- sim_landscape$V1[i]
    } 
    else {
      region_start <- sim_landscape$V2[i] 
      region_end <- sim_landscape$V2[i + 1]
      sim[region_start:region_end] <- sim_landscape$V1[i]
    }
  }
  return (as.data.frame(sim))
}

seq_len <- 30e+6
bin_sizes <- c(50e+3, 200e+3, 1e+6)

sim_land <- read.table("rho_landscape.txt")
sim_land <- expand_sim_landscape(sim_land, seq_len)
for(i in 1:length(bin_sizes)) {
  sim_land$bin <- ceiling(1:nrow(sim_land) / bin_sizes[i]) 
  binned_lands <- ddply(.data = sim_land, .variables = "bin", .fun = colMeans, .progress = "text")
  write.table(binned_lands, paste("sim.rho.", as.character(bin_sizes[i]), ".map", sep = ""), sep = "\t", quote = F)
}

sim_land <- read.table("theta_landscape.txt")
sim_land <- expand_sim_landscape(sim_land, seq_len)
for(i in 1:length(bin_sizes)) {
  sim_land$bin <- ceiling(1:nrow(sim_land) / bin_sizes[i]) 
  binned_lands <- ddply(.data = sim_land, .variables = "bin", .fun = colMeans, .progress = "text")
  write.table(binned_lands, paste("sim.theta.", as.character(bin_sizes[i]), ".map", sep = ""), sep = "\t", quote = F)
}
