# Date created: 10/11/2020
# Author: Gustavo V. Barroso
# This script writes SLiM scripts to generate the parents of our simulated trios.
# ABCDFE finishes the simulation internally, switching to viability selection

# first we create df with parameter combinations
#s = c(-1e-2, -1e-3)
#h = 0.5
#e = c(-0.2, -0.1, 0, 0.1, 0.2)
#replicate = c(1:10)

#library(tidyverse)
#df <- crossing(s, h, e, replicate)
#write.table(df, "params_df_pathways_DFE.txt", col.names = F, row.names = T, quote = F, sep = "\t")

# we will parse params_df.txt as argument for this script using bash:

# while IFS='\t' read -r -a vals; do Rscript --vanilla write_slim_script_pathways.R ${vals}; done < params_df_pathways.txt

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

idx <- as.numeric(args[1]) # which sim index
s <- as.numeric(args[2]) # which selection coefficient
h <- as.numeric(args[3]) # which dominance coefficient
e <- as.numeric(args[4]) # which epistasis coefficient
r <- as.numeric(args[5]) # which replicate

dir.create(paste("script_", idx, sep = ""))

# table generated with get_hg_data.R
exome <- read.table("~/Data/ABCDFE/human_genome/human_genome_pathways.txt", header = T)

# starts writing SLiM script
sink(paste("script_", idx, "/pathways_e_", as.character(e), "_rep_", r, "_DFE_deCODE_maps_gravel_demo.slim", sep = ""))

cat("initialize() {\n") # this is a Wright-Fisher simulation

cat(paste("setSeed(", as.character(r), ");\n", sep = "")) # for reproducibility

cat(paste("defineConstant(\"h\",", as.character(h), ");\n", sep = "")) 
#cat(paste("defineConstant(\"s\",", as.character(s), ");\n", sep = ""))
cat(paste("defineConstant(\"e\",", as.character(e), ");\n", sep = ""))

cat("initializeMutationRate(1e-8);\n")
cat("initializeMutationType(\"m0\", h, \"f\", 0);\n") # neutral	mutations
cat("initializeMutationType(\"m1\", h, \"g\", -0.01314833, 0.186);\n") # R-HSA-9675108	mutations
cat("initializeMutationType(\"m2\", h, \"g\", -0.01314833, 0.186);\n") # R-HSA-212165 mutations
cat("initializeMutationType(\"m3\", h, \"g\", -0.01314833, 0.186);\n") # R-HSA-69620 mutations
cat("initializeGenomicElementType(\"g1\", c(m1,m0), c(2.31,1.0));\n") # R-HSA-9675108	exons
cat("initializeGenomicElementType(\"g2\", c(m2,m0), c(2.31,1.0));\n") # R-HSA-212165 exons
cat("initializeGenomicElementType(\"g3\", c(m3,m0), c(2.31,1.0));\n") # R-HSA-69620 exons

cat("initializeSex(\"A\");\n\n")

# female recombination maps
cat("ratesF = c(")
for(i in 1:nrow(exome))
{
  cat(paste("asFloat(", as.character(exome$female_rate[i]), ")", sep = ""))
  if(i < nrow(exome)) { cat(",\n") }
}
cat(");\n\n")

cat("endsF = c(")
for(i in 1:nrow(exome))
{
  cat(paste("asInteger(", as.character(exome$exome_end[i]), ")", sep = ""))
  if(i < nrow(exome)) { cat(",\n") }
}
cat(");\n\n")

cat("initializeRecombinationRate(ratesF, endsF, \"F\");\n\n")

# male recombination maps
cat("ratesM = c(")
for(i in 1:nrow(exome))
{
  cat(paste("asFloat(", as.character(exome$male_rate[i]), ")", sep = ""))
  if(i < nrow(exome)) { cat(",\n") }
}
cat(");\n\n")

cat("endsM = c(")
for(i in 1:nrow(exome))
{
  cat(paste("asInteger(", as.character(exome$exome_end[i]), ")", sep = ""))
  if(i < nrow(exome)) { cat(",\n") }
}
cat(");\n\n")

cat("initializeRecombinationRate(ratesM, endsM, \"M\");\n\n")

exome_start_col <- which(names(exome) == "exome_start")
exome_end_col <- which(names(exome) == "exome_end")

for(i in seq(from = 2, to = nrow(exome), by =  2))
{
  if(exome[i, 6] == "R-HSA-9675108") {
    cat(paste("initializeGenomicElement(g1,asInteger(", as.character(exome[i, exome_start_col]), "),", "asInteger(", as.character(exome[i, exome_end_col]), "));\n", sep = ""))
  }
  else if(exome[i, 6] == "R-HSA-212165") {
    cat(paste("initializeGenomicElement(g2,asInteger(", as.character(exome[i, exome_start_col]), "),", "asInteger(", as.character(exome[i, exome_end_col]), "));\n", sep = ""))
  }
  else if(exome[i, 6] == "R-HSA-69620") {
    cat(paste("initializeGenomicElement(g3,asInteger(", as.character(exome[i, exome_start_col]), "),", "asInteger(", as.character(exome[i, exome_end_col]), "));\n", sep = ""))
  }
}
cat("}\n\n") # closing the "initialize" brackets

# Demography (Gravel model without Asia and with additional generations to increase sample (pop.) size)
cat("1 { sim.addSubpop(\"p1\", 7310); }\n\n")

cat("52080 { p1.setSubpopulationSize(14474); }\n\n")

cat("55960 {\n")
cat("\tsim.addSubpopSplit(\"p2\", 1861, p1);\n")
cat("\tp1.setMigrationRates(p2, 15e-5);\n")
cat("\tp2.setMigrationRates(p1, 15e-5);\n")
cat("}\n\n")

cat("57080 {\n")
cat("\tp2.setSubpopulationSize(1032);\n")
cat("\tp1.setMigrationRates(p2, 2.5e-5);\n")
cat("\tp2.setMigrationRates(p1, 2.5e-5);\n")
cat("}\n\n")

cat("57080:58150 {\n")
cat("\tt = sim.generation - 57080;\n")
cat("\tp2_size = round(1032 * exp(0.0038 * t));\n")
cat("\tp2.setSubpopulationSize(asInteger(p2_size));\n")
cat("}\n\n")

# Output entire parental generation (parents will undergo viability selection inside ABCDFE)
cat("58151 early() {\n")

cat("\tpotentialMothers = NULL;\n")
cat("\tpotentialFathers = NULL;\n")

cat("\tfor(indv in p2.individuals) {\n")
cat("\t\tif(indv.sex == \"F\") {\n")
cat("\t\t\tpotentialMothers = c(potentialMothers, indv);\n")
cat("\t\t}\n")
cat("\t\telse {\n")
cat("\t\t\tpotentialFathers = c(potentialFathers, indv);\n")
cat("\t\t}\n")
cat("\t}\n")
cat("\tpotentialMothers.genomes.outputVCF(\"potential_mothers.vcf\");\n")
cat("\tpotentialFathers.genomes.outputVCF(\"potential_fathers.vcf\");\n")
cat("}\n\n")

cat("58151 late() { sim.simulationFinished(); }\n\n")

# adjusting fitness effects based on genomic background (epistasis model after Desai et al 2007 Genetics)
cat("fitness(m1) {\n")
cat("n1 = genome1.countOfMutationsOfType(m1) + genome2.countOfMutationsOfType(m1);\n")
cat("if(n1 > 1) {\n")
cat("new_s = mut.selectionCoeff * (n1 - 1) ^ e;\n") # -1 for no self-interaction
cat("if(homozygous) {\n")
cat("return 1 + new_s;\n")
cat("}\n")
cat("else {\n")
cat("return 1 + h * new_s;\n")
cat("}\n")
cat("}\n")
cat("else {\n")
cat("return relFitness;\n")
cat("}\n")
cat("}\n\n")

cat("fitness(m2) {\n")
cat("n2 = genome1.countOfMutationsOfType(m2) + genome2.countOfMutationsOfType(m2);\n")
cat("if(n2 > 1) {\n")
cat("new_s = mut.selectionCoeff * (n2 - 1) ^ e;\n") # -1 for no self-interaction
cat("if(homozygous) {\n")
cat("return 1 + new_s;\n")
cat("}\n")
cat("else {\n")
cat("return 1 + h * new_s;\n")
cat("}\n")
cat("}\n")
cat("else {\n")
cat("return relFitness;\n")
cat("}\n")
cat("}\n\n")

cat("fitness(m3) {\n")
cat("n3 = genome1.countOfMutationsOfType(m3) + genome2.countOfMutationsOfType(m3);\n")
cat("if(n3 > 1) {\n")
cat("new_s = mut.selectionCoeff * (n3 - 1) ^ e;\n") # -1 for no self-interaction
cat("if(homozygous) {\n")
cat("return 1 + new_s;\n")
cat("}\n")
cat("else {\n")
cat("return 1 + h * new_s;\n")
cat("}\n")
cat("}\n")
cat("else {\n")
cat("return relFitness;\n")
cat("}\n")
cat("}\n\n")

sink()
