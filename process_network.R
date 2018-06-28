setwd('~/Desktop/CSHL/')
library(vwr)
years <- c(2012, 2013, 2014, 2015, 2016, 2017)
for(i in 1:length(years)){
  load(paste('network_', years[i], '.RData', sep = ''))
  assign(paste('network_', years[i], sep = ''), test)
  remove(test)
}

prev <- c()
for(i in 1:length(years)){
  lengths <- c()
  network <- get(paste('network_', years[i], sep = ''))
  for(j in 1:length(network)){
    lengths[j] <- length(network[[j]])
  }
  prev[i] <-  names(network[which.max(lengths)])
}

load('AA_GEN.RData')
aa_ya <- which(aa_gen[,4] == prev)
aa_seq <- c()
for(i in 1:length(prev)){
  aa_seq[i] <- which(aa_gen[,4] == prev[i])
}
prev_aa <- aa_gen[aa_seq, 5]
new_prev_aa <- c()
for(i in 1:length(prev_aa)){
  new_prev_aa[i] <- paste(unlist(prev_aa[[i]]), collapse = '')
}
hamming.neighbors(as.character(prev_aa[1]), as.character(prev_aa))
