setwd('~/Desktop/CSHL/')
load('DNA_GEN.RData')
load('AA_GEN.RData')
load('out.list_unique.years.RData')
mutations <- c('g291d', 'i156r')
years_vec <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
labels <- unlist(aa_gen[,4])
library(vwr)

mut.291 <- c()
for(i in 1:nrow(aa_gen)){
  aa <- unname(unlist(aa_gen[i,5]))
  if(length(aa) < 560){mut.291[i] <- NA & next}
  if(aa[291] == 'd'){mut.291[i] <- 'MUT'}else{mut.291[i] <- NA}
}

mut.156 <- c()
for(i in 1:nrow(aa_gen)){
  aa <- unname(unlist(aa_gen[i,5]))
  if(length(aa) < 560){mut.156[i] <- NA & next}
  if(aa[156] == 'r'){mut.156[i] <- 'MUT'}else{mut.156[i] <- NA}
}

mut.291_labels <- labels[which(mut.291 == "MUT")]
mut.156_labels <- labels[which(mut.156 == "MUT")]

mut.291_sequences <- c()
mut.156_sequences <- c()
for(i in 1:length(mut.291_labels)){
  mut.291_sequences[i] <- paste(unlist(dna_gen[which(dna_gen[,4] == mut.291_labels[i]), 5]), collapse = '')
}
for(i in 1:length(mut.156_labels)){
  mut.156_sequences[i] <- paste(unlist(dna_gen[which(dna_gen[,4] == mut.156_labels[i]), 5]), collapse = '')
}

ham_dist.291 <- function(i){
  comp <- mut.291_sequences[i]
  to <- mut.291_sequences
  ham <- hamming.neighbors(as.character(comp), as.character(mut.291_sequences))[1]
  this.name <- labels[i]
  length.ham <- length(unlist(ham))
  ham.out <- list()
  if(is.null(unlist(ham))){return(0)}else{
    for(k in 1:length.ham){
      neighbor.seq <- unlist(ham)[k]
      neighbor.which <- which(to == neighbor.seq)
      neighbor.name <- labels[neighbor.which]
      ham.out[[k]] <- neighbor.name
    }
  }
  return(unlist(ham.out))
}
ham_dist.156 <- function(i){
  comp <- mut.156_sequences[i]
  to <- mut.156_sequences
  ham <- hamming.neighbors(as.character(comp), as.character(mut.156_sequences))[1]
  this.name <- labels[i]
  length.ham <- length(unlist(ham))
  ham.out <- list()
  if(is.null(unlist(ham))){return(0)}else{
    for(k in 1:length.ham){
      neighbor.seq <- unlist(ham)[k]
      neighbor.which <- which(to == neighbor.seq)
      neighbor.name <- labels[neighbor.which]
      ham.out[[k]] <- neighbor.name
    }
  }
  return(unlist(ham.out))
}

out.list.291 <- list()
for(i in 1:length(mut.291_sequences)){
  out.list.291[[i]] <- cbind(rep(mut.291_labels[i], length(ham_dist.291(i))), ham_dist.291(i))
}
out.mat.291 <- do.call(rbind, out.list.291)
out.mat.291 <- out.mat.291[-which(out.mat.291[,2] == 0),]
out.mat.291 <- unique(out.mat.291)

out.list.156 <- list()
for(i in 1:length(mut.156_sequences)){
  out.list.156[[i]] <- cbind(rep(mut.156_labels[i], length(ham_dist.156(i))), ham_dist.156(i))
}
out.mat.156 <- do.call(rbind, out.list.156)
out.mat.156 <- out.mat.156[-which(out.mat.156[,2] == 0),]
out.mat.156 <- unique(out.mat.156)


#####PLOT 291 NETWORK
years <- c('2010', '2011', '2012', '2013', '2014')
edges <- out.mat.291
vertices <- unique(c(out.mat.291[,1], out.mat.291[,2]))
years.unique <- c()
years_colors <- c()
colors <- brewer.pal(5, "Spectral")
for(i in 1:length(vertices)){
  years.unique[i] <- unlist(dna_gen[which(dna_gen[,4] == vertices[i]),1])
  years_colors[i] <- colors[which(years == years.unique[i])]
}
frames <- years_colors 
years_colors[match(vertices, mut.291_labels)] <- "black"

graph.object <- graph.data.frame(edges, vertices, directed = F)
V(graph.object)$name <- NA
V(graph.object)$color <- years_colors
V(graph.object)$frame.color <- frames
plot(graph.object, vertex.size = 3,
     main = "291 mut network")
legend(x= -2, y = 1.002389, c(years, 'mut 291'), fill = c(colors,'black'))

#####PLOT 156 NETWORK
years <- c('2010', '2011', '2012', '2013', '2014', '2015', '2016')
edges <- out.mat.156
vertices <- unique(c(out.mat.156[,1], out.mat.156[,2]))
years.unique <- c()
years_colors <- c()
colors <- brewer.pal(7, "Spectral")
for(i in 1:length(vertices)){
  years.unique[i] <- unlist(dna_gen[which(dna_gen[,4] == vertices[i]),1])
  years_colors[i] <- colors[which(years == years.unique[i])]
}
frames <- years_colors 
years_colors[match(vertices, mut.291_labels)] <- "black"

graph.object <- graph.data.frame(edges, vertices, directed = F)
V(graph.object)$name <- NA
V(graph.object)$color <- years_colors
V(graph.object)$frame.color <- frames
plot(graph.object, vertex.size = 3,
     main = "156 mut network")
legend(x= -2, y = 1.002389, c(years, 'mut 156'), fill = c(colors,'black'))


######MAKE BARPLOT CHARTS
par(mfrow = c(2, 1))
years <- paste(aa_gen[,1], aa_gen[,2], aa_gen[3], sep = '-')
years <- as.Date(years, "%Y-%m-%d")
x.291 <- years[which(is.na(mut.291) == FALSE)]
y.291 <- rep(1, length(x.291))
x.156 <- years[which(is.na(mut.156) == FALSE)]
y.156 <- rep(1, length(x.156))
plot(x = x.291, y = y.291, pch = '|', main = "Prevalence of 291 Mutation")
plot(x = x.156, y = y.156, pch = '|', main = "Prevalence of 156 Mutation")











