args = commandArgs(TRUE)
input = as.numeric(args[1])


library(stringdist)
load('UNIQUE_DNA.RData')
dna_2015 <- dna[(which(dna[,1] == 2015 )),]
dna_2016 <- dna[(which(dna[,1] == 2016 )),]
dna_2017 <- dna[(which(dna[,1] == 2017 )),]
dna_year <- rbind(dna_2015, dna_2016, dna_2017)


dna_to <- list()
for(i in 1:dim(dna_year)[1]){
  dna_to[[i]] <- paste(unlist(dna_year[i,5]), collapse = '')
}
dna_to <- as.character(dna_to)

dist_mat <- stringdistmatrix(dna_to)

save(dist_mat, file = 'dist_mat.RData')