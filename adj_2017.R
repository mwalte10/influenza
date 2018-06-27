args = commandArgs(TRUE)
input = as.numeric(args[1])

load('UNIQUE_DNA.RData')
library(Matrix)
dna_2017 <- dna[(which(dna[,1] == 2017)),]

adj_mat_2017 <- Matrix(NA, nrow = dim(dna_2017)[1], ncol = dim(dna_2017)[1], sparse = TRUE)

for(i in 1:dim(dna_2017)[1]){
  comp <- dna_2017[,5][i]
  comp <- unlist(comp)
  for(j in 1:dim(dna_2017)[1]){
    to <- dna_2017[,5][j]
    to <- unlist(to)
    x <- length(which(comp != to))
    if(x == 1){adj_mat_2017[i,j] <- 1} 
  }
  print(i / dim(dna_2017)[1])
}

save(adj_mat_2017, file = 'adj_mat_2017.RData')