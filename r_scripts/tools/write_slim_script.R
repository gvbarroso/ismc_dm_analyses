# Date created: 10/11/2020
# Author: Gustavo V. Barroso

# table generated with get_dm_data.R
genome <- read.table("dm_tbl.txt", header = T)
genome <- genome[genome$seqnames == 1,]

# starts writing SLiM script
sink(paste("TEST_dm2L_script_1.slim", sep = ""))

cat("initialize() {\n") 

cat(paste("setSeed(1);\n", sep = "")) 

cat("defineConstant(\"N\", 20000);\n")
cat("initializeTreeSeq(simplificationRatio=100);\n")
cat("initializeMutationType(\"m1\", 0.0, \"g\", -0.01, 0.2);\n\n") # NS	mutations
cat("initializeGenomicElementType(\"g1\", m1, 1.0);\n") # exon

# comeron 2012 recombination map
#cat("ratesRec = c(")
#for(i in 1:nrow(genome))
#{
#  cat(paste("asFloat(", as.character(genome$rec_rate[i]), ")", sep = ""))
#  if(i < nrow(genome)) { cat(",\n") }
#}
#cat(");\n\n")

#cat("endsRec = c(")
#for(i in 1:nrow(genome))
#{
#  cat(paste("asInteger(", as.character(genome$end[i]), ")", sep = ""))
#  if(i < nrow(genome)) { cat(",\n") }
#}
#cat(");\n\n")

#cat("initializeRecombinationRate(ratesRec, endsRec);\n\n")
cat("initializeRecombinationRate(1e-8);\n\n")

# "fake" mutation map, just with zero in the intergenic regions to make it faster with tree-seq recording 
cat("ratesMut = c(")
for(i in 1:nrow(genome))
{
  if(!is.na(genome[i, 4])) {
    cat("asFloat(2.5e-8 * 2.31 / 3)", sep = "")
  }
  else {
    cat("asFloat(0)", sep = "")
  }
  if(i < nrow(genome)) { cat(",\n") }
}
cat(");\n\n")

cat("endsMut = c(")
for(i in 1:nrow(genome))
{
  cat(paste("asInteger(", as.character(genome$end[i]), ")", sep = ""))
  if(i < nrow(genome)) { cat(",\n") }
}
cat(");\n\n")

cat("initializeMutationRate(ratesMut, endsMut);\n\n")

start_col <- which(names(genome) == "start")
end_col <- which(names(genome) == "end")

for(i in 1:nrow(genome))
{
  cat(paste("initializeGenomicElement(g1,asInteger(", as.character(genome[i, start_col]), "),",
            "asInteger(", as.character(genome[i, end_col]), "));\n", sep = ""))
}
cat("}\n\n") # closing the "initialize" brackets

cat("s1 20000 late() {\n")
cat("sim.treeSeqOutput(\"dm2L_bgs.trees\");\n")
cat("sim.simulationFinished();\n")
cat("}\n\n")

cat("1 {\n")
cat("sim.addSubpop(\"p1\", N);\n")
cat("sim.rescheduleScriptBlock(s1, start=10*N, end=10*N);\n")
cat("}\n\n")

sink()
