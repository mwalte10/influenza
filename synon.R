setwd('~/Desktop/CSHL/')

# dna <- read.fasta('DNA.fasta')
# aa <- read.fasta('AA.fasta')
load('id.RData')
load('DNA.RData')
load('AA.RData')

unique_dna <- unique(dna_temp)
unique_aa <- unique(aa_temp)

##this isn't working but need to do something like this
counts_dna <- c()
for(i in 1:length(unique_dna)){
  counts_dna[i] <- length(which(dna_temp == unique_dna[i]))
}