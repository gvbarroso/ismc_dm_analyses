# Created: 07/08/2017
# Last modified: 07/08/2018
# Author: Gustavo Barroso
# This script performs simulations of recombination and mutation landscapes under diverse scenarios

##############################################
#
# General Parameters 
#
##############################################

sequence_length <- 30e+6
N0 <- 1e+5
Ne <- N0 # used when simulating fluctuating pop. sizes
little_r <- 1e-9 # recombination rate per site per generation
mean_rho <- 4 * N0 * little_r
little_mu <- 2e-8 # mutation rate per site per generation
mean_theta <- 4 * N0 * little_mu
num_haploids <- 10 # sample size (in haploids)


##############################################
#
# Demographic History 
#
##############################################

# For demographic history, the directory stays the same

# Introgression (1-pulse event between 2 populations)
admix <- FALSE
if(admix) {
  admix_proportion <- 0.1 # proportion of genetic material from source to target population
  admix_time <- 0.125 # time (in coalescent units) of secondary contact between populations
  split_time <- 2.0 # time (in coalescent units) of separation between populations
}  

# discretized time intervals used by iSMC
number_of_intervals <- 30
time_boundaries <- rep(0, number_of_intervals)
for(i in 2:number_of_intervals) {
  time_boundaries[i] = - log(1 - ((i - 1) / number_of_intervals))
}

# vector of time points (in coal. units) where demographic changes occur
time_demo_changes <- c(0.18232156, 0.83624802)
# by how much the pop. size changes (past-wards at each step) at the above time points 
fold_changes <- c(1, 0.285714285714286, 1.5)

# Now, since we have fluctuating pop. sizes, we must ajust N0 to have the (time-average) Ne that we want
# computes total (maximum) number of generations by assuming max. scaled coal. time is equal to max. time in discretisation
max_time <- time_boundaries[number_of_intervals] # can be viewed as ~expected maximum tree height in the ARG
num_gen_total <- 4 * Ne * max_time

# number of generations spent in-between each pop. size change
num_gen_ib <- 4 * Ne * time_demo_changes
num_gen_vec <- num_gen_ib[1]
if(length(num_gen_ib) > 1) {
  for(i in 2:length(num_gen_ib)) {
    num_gen_vec[i] <- num_gen_ib[i] - num_gen_ib[i - 1]
  }
}
num_gen_vec <- c(num_gen_vec, num_gen_total - num_gen_vec[length(num_gen_vec)])

# we solve for N0 using the harmonic mean formula:
harmonic_pop_sizes <- 0
for(i in 1:length(num_gen_vec)) {
  harmonic_pop_sizes <- harmonic_pop_sizes + num_gen_vec[i] / fold_changes[i]
}
N0 <- (harmonic_pop_sizes * Ne) / num_gen_total
  

# since SCRM sets theta and rho based on N0, we must update these values
# IMPORTANT: such rates are "genome-wide 'means'" (relative to mu and little_r), but NOT "time-wise 'means'" (relative to Ne). 
mean_theta <- 4 * N0 * little_mu 
mean_rho <- 4 * N0 * little_r

# plots demo. history
coal_times <- c(0.0, time_demo_changes, time_boundaries[number_of_intervals])
generation_times <- 4 * N0 * coal_times
pop_sizes <- c(fold_changes, fold_changes[length(fold_changes)])
pop_sizes <- pop_sizes * N0

pdf("demographic_history.pdf")
plot(x = generation_times, 
     y = pop_sizes, type = "s", lwd = 2.5, xaxt = "n", ylab = "Ne",  xlab = "Time (4Ne Generations)")
# plot x-axis to match time boundaries (best if x-axis is plotted in log-scale):
x_axis <- round(4 * N0 * time_boundaries)
# plot x axis in equally spaced boundaries (best if x-axis is plotted in linear scale)
x_axis <- round(4 * N0 * seq(from = 0.0, to = max(time_boundaries), length.out = 40))
axis(1, at = x_axis, las = 2, cex.axis = 0.5)
dev.off()


##############################################
#
# Gamma model of spatial variation in Rho
#
##############################################

# Gamma distribution of rho
alpha_rho <- 0.5
beta_rho <- alpha_rho

# frequency of change in rho values along the genome 
rho_transition_prob <- 5e-5

number_of_rho_transitions = 0
rho_transition_points = c(NULL)
current_transition_point = 0
rho_span_vector = c(NULL)

# gets the points where rho values change
while(current_transition_point < sequence_length) {
  number_of_rho_transitions <- number_of_rho_transitions + 1
  rho_span <- rgeom(1, rho_transition_prob)
  rho_span_vector <- c(rho_span_vector, rho_span)
  current_transition_point <- current_transition_point + rho_span
  rho_transition_points[number_of_rho_transitions] <- current_transition_point
}

# if the last transition point happens to fall outside the range of our sequence, we delete it
if(rho_transition_points[number_of_rho_transitions] > sequence_length){
  rho_transition_points <- rho_transition_points[- number_of_rho_transitions]
  number_of_rho_transitions <- number_of_rho_transitions - 1
  rho_span_vector <- rho_span_vector[- number_of_rho_transitions]
}

# gets the rho values from the gamma distribution
rho_values <- rgamma(n = number_of_rho_transitions + 1, shape = alpha_rho, rate = beta_rho)
rho_values <- append(rho_values, rho_values[length(rho_values)], after = length(rho_values))
rho_values <- rho_values * mean_rho # scale by genome-wide average rho to get a meaningful rate
first_rho <- rho_values[1]

# prints and writes rho landscape to file
rho_transition_points <- append(rho_transition_points, 0, after = 0)
rho_transition_points <- append(rho_transition_points, sequence_length, after = length(rho_transition_points))
# writes landscapes scaled back to (time-average) Ne, since this is how iSMC infers it
rho_landscape <- cbind(as.data.frame(rho_values * (Ne / N0)), as.data.frame(rho_transition_points))

pdf("sim_rho_landscape.pdf")
plot(y = rho_landscape$rho_values, x = as.integer(rho_landscape$rho_transition_points / 1e+3),
     type = "s", ylab = expression(rho), xlab = "Position (kb)", lwd = 0.5, main = "Gamma-simulated recombination landscape")
abline(h = mean_rho * (Ne / N0), lty = 2, col = "red")
dev.off()
write.table(rho_landscape, file = "rho_landscape.txt", quote = F, row.names = F, col.names = F)

# writes simulated parameters to file 
# we write mean_theta and mean_rho as a function of Ne, not N0, so we re-scale them
values <- c(mean_theta * (Ne / N0), mean_rho * (Ne / N0), alpha_rho, (1 - rho_transition_prob), little_r, N0, Ne)
names <- c("mean_theta", "mean_rho", "rho.alpha", "r_ii", "little_r", "N0", "Ne") 
sim_params <- cbind(names, values)
write.table(sim_params, file = "sim_params.txt", quote = F, row.names = F, col.names = F, sep = "\t")

##############################################
#
# Gamma model of spatial variation in Theta
#
##############################################

# Gamma distribution of theta values 
alpha_theta <- 2.5
beta_theta <- alpha_theta

# frequency of change in theta along the sequence
theta_transition_prob <- 2e-5

number_of_theta_transitions = 0
theta_transition_points = c(NULL)
current_transition_point = 0
theta_span_vector = c(NULL)

# gets the points where theta values change
while(current_transition_point < sequence_length) {
  number_of_theta_transitions <- number_of_theta_transitions + 1
  theta_span <- rgeom(1, theta_transition_prob)
  theta_span_vector <- c(theta_span_vector, theta_span)
  current_transition_point <- current_transition_point + theta_span
  theta_transition_points[number_of_theta_transitions] <- current_transition_point
}

# if last transition != sequence length, we delete it
if(theta_transition_points[number_of_theta_transitions] >= sequence_length) {
  theta_transition_points <- theta_transition_points[- number_of_theta_transitions]
  number_of_theta_transitions <- number_of_theta_transitions - 1
  theta_span_vector <- theta_span_vector[- number_of_theta_transitions]
}

theta_values <- numeric(length = number_of_theta_transitions + 1)

theta_values <- rgamma(n = number_of_theta_transitions + 1, shape = alpha_theta, rate = beta_theta)

theta_values <- append(theta_values, theta_values[length(theta_values)], after = length(theta_values))
# transforms to mean 1
theta_values <- theta_values * mean_theta
first_theta <- theta_values[1]

# assembles the theta landscape for reference
theta_transition_points <- append(theta_transition_points, 0, after = 0)
theta_transition_points <- append(theta_transition_points, sequence_length, after = length(theta_transition_points))
theta_landscape <- as.data.frame(cbind(theta_values * (Ne / N0),theta_transition_points))

# go back
pdf("sim_theta_landscape.pdf")
plot(y = as.numeric(theta_landscape[,1]), x = as.numeric(theta_landscape[,2]), lwd = 2,
     type = "s", ylab = expression(theta), xlab = "Position (bp)", main = "Theta landscape")
abline(h = mean_theta, lty = 2, col = "blue")
dev.off()

write.table(theta_landscape, file = "theta_landscape.txt", quote = F, row.names = F, col.names = F)

# writes simulated parameters to file
values <- c(alpha_theta, theta_transition_prob, little_mu, N0)
names <- c("theta.alpha", "t_ij", "little_mu", "N0")
sim_params <- cbind(names, values)
write.table(sim_params, file = "sim_params_theta.txt", quote = F, row.names = F, col.names = F, sep = "\t")


##############################################
#
# Creating the SCRM command line
#
##############################################

# build an SCRM commandline based on the parameters used so far in the script
sink("scrm_commandline.txt")
if(exists("rho_landscape")) {
  cat(paste(as.character(num_haploids), "1 -r", as.character(first_rho * sequence_length),
            as.character(sequence_length), "", sep = " "))
  for (i in 2:(number_of_rho_transitions + 1)){
    cat("-sr ")
    cat(as.character(rho_transition_points[i]))
    cat(" ")
    cat(as.character(rho_values[i] * sequence_length))
    cat(" ")
  }
} else {
  cat(paste(as.character(num_haploids), "1 -r", as.character(mean_rho * sequence_length),
            as.character(sequence_length), "", sep = " "))
}
if(exists("theta_landscape")) {
  cat(paste("-t", as.character(first_theta * sequence_length), "", sep = " "))
  for(i in 2:(number_of_theta_transitions + 1)) {
    cat("-st ")
    cat(as.character(theta_transition_points[i]))
    cat(" ")
    cat(as.character(theta_values[i] * sequence_length))
    cat(" ")
  }
} else {
  cat(paste("-t", as.character(mean_theta * sequence_length)))
  cat(" ")
}
if(exists("time_demo_changes")) {
  for(i in 1:length(time_demo_changes)){
    cat("-eN ")
    cat(as.character(time_demo_changes[i]))
    cat(" ")
    cat(as.character(fold_changes[i]))
    cat(" ")
  }
}
if(exists("admix_proportion")) {
  cat("-es ")
  cat(as.character(admix_time))
  cat(" ")
  cat("1")
  cat(" ")
  cat(as.character(1 - admix_proportion))
  cat(" ")
  cat("-ej ")
  cat(as.character(split_time))
  cat(" ")
  cat("1")
  cat(" ")
  cat("2")
}
cat("-T") 
cat("\n")
sink()

##############################################
#
# Running SCRM 
#
##############################################

# reads generated SCRM commandline to map mutations into haplotypes
require(scrm)
require(stringr)

# This part of the script is compatible with R version > 3.6.0

command_line <- readLines("scrm_commandline.txt")

# number of ARGs to be generated (replicates) based on the simulated landscape
num_reps <- 10
# creates one dir per replicate
for(i in 1:num_reps) {
  rep_dir <- paste("rep_", as.character(i), sep = "")
  if(!dir.exists(rep_dir)){
    dir.create(rep_dir, showWarnings = F)
  } else {
    print(paste("Dir", rep_dir, "already exists! Over-writing..."))
    dir.create(rep_dir, showWarnings = F)
  }
}

pb <- txtProgressBar(min = 0, max = num_reps, style = 3)
for(i in 1:num_reps) {
  setTxtProgressBar(pb, i)
  set.seed(i)
  
  command_line <- readLines("scrm_commandline.txt")
  
  full_arg <- scrm(command_line) # runs scrm
  
  lapply(full_arg$trees[1], write, paste("rep_", i, "/rep_", i, ".newick", sep = ""), append = T)
  
  haps <- as.data.frame(t(as.data.frame(full_arg$seg_sites))) # organises haplotypes
  pos <- row.names(haps) 
  pos <- str_replace(pos, "X", "")
  # gets infinite-sites SNPs positions
  num_pos <- as.numeric(str_replace(pos, "e.", "e-")) 
  # converts to finite loci coordinates
  snp_loci <- round(num_pos *  sequence_length, 0)
  # moves multiple hits to the right
  snp_loci[which(duplicated(snp_loci))] <- snp_loci[which(duplicated(snp_loci))] + 1
  # if it happ  ens again, we discard them
  if(length(which(duplicated(snp_loci))) > 0) {
    dup_pos <- which(duplicated(snp_loci))
    snp_loci <- snp_loci[-dup_pos]
    haps <- haps[-dup_pos,]
  }
  haps$pos <- snp_loci
  
  full_haplotyptes <- as.data.frame(matrix(ncol = num_haploids, nrow = sequence_length))
  full_haplotyptes[is.na(full_haplotyptes)] <- 0
  # maps mutations into haplotypes
  for(j in 1:num_haploids) {
    names(full_haplotyptes)[j] <- paste("hap_", as.character(j), sep = "")
    hap_seq <- numeric(length = sequence_length)
    # gets snp coordinates for haplotype j
    snp_pos <- haps[which(haps[,j] == 1), ncol(haps)]
    hap_seq[snp_pos] <- 1
    full_haplotyptes[,j] <- hap_seq
  }
  file_path <- paste(getwd(), "/rep_", as.character(i), "/rep_", as.character(i), sep = "")
  # writes with col. names as a header
  write.table(full_haplotyptes, file = file_path, row.names = F, col.names = T)
}
close(pb)

