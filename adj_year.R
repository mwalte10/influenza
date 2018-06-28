args = commandArgs(TRUE)
input = as.numeric(args[1])
years <- c(2012, 2013, 2014, 2015, 2016, 2017)
years <- years[input]

#setwd('~/Desktop/CSHL/')
load('UNIQUE_DNA.RData')
library(Matrix)
dna_year <- dna[(which(dna[,1] == years)),]

#assign(paste0('adj_mat_', years, sep = ''), matrix(0, nrow = dim(dna_year)[1], ncol = dim(dna_year)[1])) 
mat <- matrix(0, nrow = dim(dna_year)[1], ncol = dim(dna_year)[1])

for(i in 1:dim(dna_year)[1]){
  comp <- dna_year[,5][i]
  comp <- unlist(comp)
  for(j in 1:dim(dna_year)[1]){
    to <- dna_year[,5][j]
    to <- unlist(to)
    x <- length(which(comp != to))
    ifelse(x == 1, mat[i,j] <-  1, 
           mat[i,j] <-  0)
  }
  print(i / dim(dna_year)[1])
}

colnames(mat) <- unlist(dna_year[,4])
rownames(mat) <- unlist(dna_year[,4])


counts <- colSums(mat)
names(counts) <- unlist(dna_year[,4])

###Essentially gets the coordinates of the ones which are interacting
cols <- (which(mat==1, arr.ind=TRUE))[,2]
rows <- (which(mat==1, arr.ind=TRUE))[,1]
new <- cbind(unlist(dna_year[,4])[cols], unlist(dna_year[,4])[rows])
test <- network(new, matrix.type = 'edgelist', directed = FALSE)
plot.network(test, main = "2013 DNA Sequence Network \nAll Unique Sequences", vertex.col = 'grey', 
             usecurve = 0, #vertex.cex = seq(0.2, 1.3, length.out = max(counts))
             jitter= FALSE)

sparse_mat <- as(mat, "sparseMatrix") 
save(test, file = paste('network_', years, '.RData'))
save(sparse_mat, file = paste('adj_mat_', years, '.RData'))
