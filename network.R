library(network)
library(vrw)

load('UNIQUE_DNA.RData')
years_vec <- c(2012, 2013, 2014, 2015, 2016, 2017)
par(mfrow = c(2,3))

years <- years_vec[input]
dna_year <- dna[(which(dna[,1] == years)),]


dna_to <- list()
for(i in 1:dim(dna_year)[1]){
  dna_to[[i]] <- paste(unlist(dna_year[i,5]), collapse = '')
}
dna_to <- as.character(dna_to)
dna_year_id <- unlist(dna_year[,4])

ham_dist_one <- function(i){
  comp <- dna_to[i]
  to <- dna_to
  ham <- hamming.neighbors(as.character(comp), as.character(dna_to))[1]
  if(is.null(unlist(ham))){return(0)}else{
    ham.out <- c()
    for(j in 1:length(unlist(ham))){
      ham.out[j] <-  which(dna_to == unlist(ham)[j])
    }
    return(ham.out)
  }
}

test <- list()
for(i in 1:length(dna_to)){
  test[[i]] <- ham_dist_one(i)
  print(i / length(dna_to))
}
names(test) <- dna_year_id

left <- list()
right <- list()
for(i in 1:length(dna_to)){
  # if((unlist(test[i]) != 0)[1]){left[[i]] <- rep(names(test[i]), length(test[[i]]))}
  # if((unlist(test[i]) != 0)[1]){right[[i]] <- dna_year_id[unlist(unname(test[i]))]}
  if(length((unlist(test[i]) != 0)) > 2){left[[i]] <- rep(names(test[i]), length(test[[i]]))}
  if(length((unlist(test[i]) != 0)) > 2){right[[i]] <- dna_year_id[unlist(unname(test[i]))]}
}
left <- unlist(left)
right <- unlist(right)
matrix_network <- cbind(left, right)
mat.net <- network(matrix_network, matrix.type = 'edgelist', directed = FALSE)
plot.network(mat.net, main = paste(years_vec[input], 'DNA Sequence Network \nMore Than One Neighbor', sep = ' '), vertex.col = 'grey', 
             usecurve = 0, #vertex.cex = seq(0.2, 1.3, length.out = max(counts))
             jitter= FALSE)

