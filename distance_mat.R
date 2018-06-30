args = commandArgs(TRUE)
input = as.numeric(args[1])

library(network)
library(vwr)
library(RColorBrewer)
load('UNIQUE_DNA.RData')


dna_to <- list()
for(i in 1:dim(dna_year)[1]){
  dna_to[[i]] <- paste(unlist(dna_year[i,5]), collapse = '')
}
dna_to <- as.character(dna_to)

stringdistmatrix(dna_to)