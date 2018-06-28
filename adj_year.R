args = commandArgs(TRUE)
input = as.numeric(args[1])
years <- c(2012, 2013, 2014, 2015, 2016, 2017)
years <- years[input]

load('UNIQUE_DNA.RData')
library(Matrix)
dna_year <- dna[(which(dna[,1] == years)),]

paste('adj_mat_', years, sep = '') <- Matrix(NA, nrow = dim(dna_year)[1], ncol = dim(dna_year)[1], sparse = TRUE)

for(i in 1:dim(dna_year)[1]){
  comp <- dna_year[,5][i]
  comp <- unlist(comp)
  for(j in 1:dim(dna_year)[1]){
    to <- dna_year[,5][j]
    to <- unlist(to)
    x <- length(which(comp != to))
    if(x == 1){adj_mat_year[i,j] <- 1} 
  }
  print(i / dim(dna_year)[1])
}

save(paste('adj_mat_', years, sep = ''), 
     file = paste('adj_mat_', years, '.RData', sep = ''))