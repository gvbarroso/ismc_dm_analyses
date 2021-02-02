# Date created: 02/12/2019
# Author: Gustavo V. Barroso
# Retrieves and organizes data from droso Genome, to set up SLiM simulations

library("biomaRt")
library("GenomicRanges")
library("tidyverse")
library("plyr")
library("magrittr")
library("BRGenomics")

setwd("~/Data/iSMC/theta_paper/slim_sims/")

#######################
#
# Loads & Filters data
#
#######################

chr_labels <- c("2L", "2R", "3L", "3R")

# gets exome data from drosophila genome (accessed on January 25 2021 to perform simulations)
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl") 

droso_exome_ensembl <- getBM(mart = ensembl, filters = "chromosome_name", values = chr_labels,
                             attributes = c("ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end"))

droso_exome <- dplyr::arrange(.data = droso_exome_ensembl, chromosome_name, exon_chrom_start)

# filters out overlapping exons / genes
droso_exome <- ddply(.data = droso_exome, .variables = "chromosome_name", .fun = distinct, exon_chrom_start, .keep_all = T)
droso_exome <- ddply(.data = droso_exome, .variables = "chromosome_name", .fun = distinct, exon_chrom_end, .keep_all = T)

droso_exome_gr <- makeGRangesFromDataFrame(droso_exome, keep.extra.columns = T)

# iteratively gets rid of the final overlapping exons
curr_num_ranges <- length(ranges(droso_exome_gr))
prev_num_ranges <- curr_num_ranges + 1
while(prev_num_ranges > curr_num_ranges) {
  droso_exome_gr <- droso_exome_gr[unique(findOverlaps(droso_exome_gr, type = "any", select = "first"))];
  prev_num_ranges <- curr_num_ranges
  curr_num_ranges <- length(ranges(droso_exome_gr))
}

# merges contiguous exons from different genes 
cont_gr <- as.data.frame(droso_exome_gr)
good_rows <- rep(T, nrow(cont_gr))
for(i in 2:nrow(cont_gr)) {
  if(cont_gr[i, 2] == (cont_gr[i - 1, 3] + 1)) {
    cont_gr[i - 1, 3] <- cont_gr[i, 3]
    good_rows[i] <- F
  }
}
cont_gr <- cont_gr[good_rows, c(1:3, 6)]
row.names(cont_gr) <- 1:nrow(cont_gr)

cont_gr$seqnames <- as.numeric(cont_gr$seqnames)
cont_gr$start <- as.numeric(cont_gr$start)
cont_gr$end <- as.numeric(cont_gr$end)

# updates GR
droso_exome_gr <- makeGRangesFromDataFrame(cont_gr, keep.extra.columns = T)

# rec map
rec_map <- read.table("rec_map_comeron.bedgraph", header = T, sep = "\t")
names(rec_map) <-  c("chrom", "chromStart", "chromEnd", "rate")

rec_map$chromStart <- as.numeric(rec_map$chromStart)
rec_map$chromEnd <- as.numeric(rec_map$chromEnd)
rec_map$rate <- as.numeric(rec_map$rate)

rec_map <- dplyr::arrange(.data = rec_map, chrom, chromStart)
rec_map$rate <- rec_map$rate * 1e-8 # converts from cM/Mb to r
rec_map_gr <- makeGRangesFromDataFrame(rec_map, keep.extra.columns = T)

# Organizes exonic / non-exonic data
droso_exome_df <- as.data.frame(droso_exome_gr)

# complement of the exome
non_exon_gr <- GenomicRanges::gaps(droso_exome_gr)
non_exon_df <-  as.data.frame(non_exon_gr)
non_exon_df$ensembl_exon_id <- NA

# merges exonic and non-exonic regions
dm_genome_gr <- rbind.data.frame(droso_exome_df, non_exon_df)  %>%
  dplyr::arrange(seqnames, start) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

hits <- findOverlaps(query = dm_genome_gr, subject = rec_map_gr) 

dm_genome_df <- as.data.frame(dm_genome_gr)
dm_genome_df <- dm_genome_df[-(4:5)]
dm_genome_df$rec_window <- NA
dm_genome_df[queryHits(hits), 5] <- subjectHits(hits)
dm_genome_df$rec_rate <- 0

for(i in 1:nrow(dm_genome_df)) {
  if(!is.na(dm_genome_df[i, 5])) {
    dm_genome_df[i, 6] <- rec_map[dm_genome_df[i, 5], 4]
  }
}

# accumulates coordinates bc slims treats everything as 1 chr
chr_ends <- dm_genome_df[which(dm_genome_df$seqnames != dplyr::lag(dm_genome_df$seqnames)) - 1, 3]

dm_genome_df[dm_genome_df$seqnames == 2, 2] <- dm_genome_df[dm_genome_df$seqnames == 2, 2] + chr_ends[1]
dm_genome_df[dm_genome_df$seqnames == 2, 3] <- dm_genome_df[dm_genome_df$seqnames == 2, 3] + chr_ends[1]
dm_genome_df[dm_genome_df$seqnames == 3, 2] <- dm_genome_df[dm_genome_df$seqnames == 3, 2] + chr_ends[1] + chr_ends[2]
dm_genome_df[dm_genome_df$seqnames == 3, 3] <- dm_genome_df[dm_genome_df$seqnames == 3, 3] + chr_ends[1] + chr_ends[2]
dm_genome_df[dm_genome_df$seqnames == 4, 2] <- dm_genome_df[dm_genome_df$seqnames == 4, 2] + chr_ends[1] + chr_ends[2] + chr_ends[3]
dm_genome_df[dm_genome_df$seqnames == 4, 3] <- dm_genome_df[dm_genome_df$seqnames == 4, 3] + chr_ends[1] + chr_ends[2] + chr_ends[3]

# updates
chr_ends <- dm_genome_df[which(dm_genome_df$seqnames != dplyr::lag(dm_genome_df$seqnames)) - 1, 3]

# cleaning
dm_genome_df <- dm_genome_df[-5]
dm_genome_df$seqnames <- as.numeric(dm_genome_df$seqnames)

# adds free rec. at chr boundaries 
dm_genome_df <- dm_genome_df %>% add_row(seqnames = 1, start = chr_ends[1] + 1, end = chr_ends[1] + 2, rec_rate = 0.5,
                                         .before = which(dm_genome_df$seqnames != dplyr::lag(dm_genome_df$seqnames))[1])

dm_genome_df <- dm_genome_df %>% add_row(seqnames = 2, start = chr_ends[2] + 1, end = chr_ends[2] + 2, rec_rate = 0.5,
                                         .before = which(dm_genome_df$seqnames != dplyr::lag(dm_genome_df$seqnames))[2])

dm_genome_df <- dm_genome_df %>% add_row(seqnames = 3, start = chr_ends[3] + 1, end = chr_ends[3] + 2, rec_rate = 0.5,
                                         .before = which(dm_genome_df$seqnames != dplyr::lag(dm_genome_df$seqnames))[3])

write.table(dm_genome_df, "~/Data/iSMC/theta_paper/slim_sims/dm_tbl.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# creates a table that will be used to format the chr column of SLiM-output VCF files 
chr_labels

chr_ends <- c(dm_genome_df[which(dm_genome_df$seqnames != dplyr::lag(dm_genome_df$seqnames)) - 1, 3], max(dm_genome_df$end))

chr_limits <- cbind.data.frame(chr_labels, ends_list) # will be used by format_vcf.R || filter_vcfs.sh
write.table(chr_limits, "~/Data/iSMC/theta_paper/slim_sims/chr_limits.txt", row.names = F, col.names = T, quote = F, sep = "\t")

