setwd('~/Desktop/CSHL/')
load('DNA_GEN.RData')
load('out.list_unique.years.RData')
#load('out.list_unique.RData')
#load('AA_GEN.RData')
library(vwr)
library(igraph)
library(RColorBrewer)
years_vec <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
load("preload.RData")
dna_to <- important_list[[1]]
labels <- important_list[[2]]
years.labels <- important_list[[3]]
remove(important_list)

lengths <- c()
for(i in 1:length(out.list)){
  lengths[i] <- dim(out.list[[i]])[1]
}

top_five <- tail(sort(lengths), 5)
network_top_five <- list()
for(i in 1:length(top_five)){
  network_top_five[[i]] <- out.list[[which(top_five[i] == lengths)]]
}

edges <- do.call(rbind, network_top_five)
edges <- edges[,1:2]


expanded.list <- list()
for(i in 1:length(network_top_five)){
  net <- network_top_five[[i]]
  neighbor.list <- list()
  for(j in 1:dim(net)[1]){
    if(length(which(labels == net[j,2])) == 0){next}
    neighbor.list[[j]] <- out.list[[which(labels == unname(net[j,2]))]]
  }
  expanded.list[[i]] <- do.call(rbind, neighbor.list)
}
expanded.network <- do.call(rbind, expanded.list)
expanded.edges <- expanded.network[,1:2]
edges.new <- rbind(edges, expanded.edges)
vertices.new <- unique(c(edges.new[,1], edges.new[,2]))

years <- c()
for(i in 1:length(vertices.new)){
  years[i] <- dna_gen[which(dna_gen[,4] == vertices.new[i]), 1]
  print(i / length(vertices.new))
}
years <- unlist(years)
years.colors <- c()
colors <- c(brewer.pal(8, "Accent"), "red")
for(i in 1:length(years)){
  years.colors[i] <- colors[which(years_vec == years[i])]
}
graph.object <- graph.data.frame(edges.new, vertices.new, directed = F)
V(graph.object)$color <- years.colors
plot(graph.object, vertex.size = 3, vertex.label =NA)



#####FIND CONNECTIONS
labels_top_five <- c()
for(i in 1:length(top_five)){
  labels_top_five[i] <- labels[which(top_five[i] == lengths)]
}

seqs_top_five <- c()
for(i in 1:length(top_five)){
  seqs_top_five[i] <- dna_to[which(top_five[i] == lengths)]
}
  
ham_dist <- function(j){
  comp <- seqs_top_five[j]
  to <- seqs_top_five
  ham <- hamming.neighbors(as.character(comp), as.character(seqs_top_five))
  dist <- list()
  for(i in 1:length(ham)){
    dist[[i]] <- c(labels_top_five[which(ham[[i]] == seqs_top_five)], i)
  }
  return(dist)
}
  

extend_edges_from <- c(labels_top_five[1], rep(labels_top_five[2], 2), rep(labels_top_five[3], 2), labels_top_five[4], rep(labels_top_five[5], 2)) 
extend_edges_to <- c(labels_top_five[4], labels_top_five[5], labels_top_five[3], labels_top_five[5],
                     labels_top_five[2], labels_top_five[1], labels_top_five[2], labels_top_five[3])
distances <- c(2, 1, 10, 9, 10, 2, 1, 9)
connections <- cbind(extend_edges_from, extend_edges_to, distances)
connections <- connections[5:8,]
  
########
edges.new <- cbind(edges.new, rep(1, nrow(edges.new)))
edges.new <- rbind(edges.new, connections)
vertices.new <- unique(c(edges.new[,1], edges.new[,2]))

years <- c()
for(i in 1:length(vertices.new)){
  years[i] <- dna_gen[which(dna_gen[,4] == vertices.new[i]), 1]
  print(i / length(vertices.new))
}
years <- unlist(years)
years.colors <- c()
colors <- c(brewer.pal(8, "Accent"), "red")
for(i in 1:length(years)){
  years.colors[i] <- colors[which(years_vec == years[i])]
}
graph.object <- graph.data.frame(edges.new, vertices.new, directed = F)
V(graph.object)$color <- years.colors
#E(graph.object)$width <- edges.new[,3]
plot(graph.object, vertex.size = 3, vertex.label =NA)
legend(-1.514017, 0.761, years_vec, fill = colors)
