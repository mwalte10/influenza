###SET UP INFO
{
setwd('~/Desktop/CSHL/data')
load('DNA_GEN.RData')
library(vwr)
library(igraph)
library(RColorBrewer)
library(DiagrammeR)
years_vec <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
load("preload.Rdata")
dna_to <- important_list[[1]]
labels <- important_list[[2]]
years.labels <- important_list[[3]]
remove(important_list)
}

##IF NEEDED TO REBUILD MAT
{
#####BUILD ADJ MAT
  
dna_to <- c()
for(i in 1:length(labels)){
  seq <- unlist(dna_gen[which(dna_gen[,4] == labels[i]),5])[1: (566 * 3)]
  if(any(is.na(seq))){next}
  dna_to[i] <- paste(seq, collapse = '')
  print(i / length(labels))
}
  
  ham_dist <- function(j){
    comp <- dna_to[j]
    to <- dna_to
    ham <- unlist(hamming.neighbors(as.character(comp), to)[1])
    if(is.null(unlist(ham))){return(0)}
    outputs <- c()
    for(i in 1:length(ham)){
      outputs[i] <- which(dna_to == ham[i])
    }
    return(outputs)
  }
  
  mat <- matrix(rep(0, length(labels)^2), ncol = length(labels), nrow = length(labels))
  for(j in 1:length(labels)){
    mat[j,ham_dist(j)] <- 1
    print(j / length(labels))
  }
#remove those without any neighbors
remove <- which(colSums(mat) == 0)
labels.new <- labels[-remove]
dna_to.new <- dna_to[-remove]
mat.new <- mat[-remove, -remove]

###remove all the numbers that end up in neighbors
one_neighbor <- which(colSums(mat.new) == 1)
neighbors <- c()
for(i in 1:length(one_neighbor)){
  x <- which(mat.new[,one_neighbor[i]] == 1)
  if(length(which(one_neighbor == x)) == 0){neighbors[i] <- NA}else{neighbors[i] <- one_neighbor[i]}
}

remove <- neighbors[which(is.na(neighbors) == FALSE)]
labels.new.new <- labels.new[-remove]
dna_to.new.new <- dna_to.new[remove]
mat.new.new <- mat.new[-remove, -remove]


save(labels.new.new, file = 'labels_for_simp.adj.mat.RData')
save(mat.new.new, file = 'adj_mat.RData')
#####NEED TO FIND NEW LABELS

#######USING IGRAPH
graph.mat <- graph.adjacency(mat.new.new, mode = "undirected")
save(graph.mat, file = "graph.mat.RData")

###find colors for vertices
{years <- c()
for(i in 1:length(labels)){
  years[i] <- dna_gen[which(dna_gen[,4] == labels[i]), 1]
  print(i / length(labels))
}
years <- unlist(years)
years.colors <- c()
colors <- c(brewer.pal(8, "Accent"), "red")
for(i in 1:length(years)){
  years.colors[i] <- colors[which(years_vec == years[i])]
}
graph.mat <- graph.adjacency(mat.new.new, mode = "undirected")
V(graph.mat)$name <- NA
V(graph.mat)$color <- years.colors
plot(graph.mat, vertex.name = NA, vertex.size = 2, 
     main = "Network Data, Top Down", layout = layout_with_kk(graph.mat))
legend(-1.662914, 0.85, years_vec, fill = colors)}
}

##START FROM HERE IF MAT IS OK
{
load('labels_for_simp.adj.mat.RData')
load('graph.mat.RData')
load('dfs.analysis.RData')
load('colors_for_igraph.RData')
load('adj_mat.RData')

V(graph.mat)$name <- NA
V(graph.mat)$color <- years.colors
plot(graph.mat, vertex.name = NA, vertex.size = 2, 
     main = "Network Data, Top Down", layout = layout_nicely(graph.mat))
colors <- c(brewer.pal(8, "Accent"), "red")
legend(-1.662914, 0.85, years_vec, fill = colors)
}

##FIND THE LARGEST (n.cluster) AMOUNT OF CLUSTERS
{
clusters <- clusters(graph.mat)
n.cluster <- 1
biggest_clusters <- tail(sort(clusters$csize), n.cluster)
group <- list()
for(i in 1:length(biggest_clusters)){
  group[[i]] <- which(clusters$csize == biggest_clusters[i])
}
group <- unlist(group)
clusters_groups <- list()
for(i in 1:length(group)){
  clusters_groups[[i]] <- which(clusters$membership == group[i])
}
labels_vec <- labels.new.new[unlist(clusters_groups)]
save(labels_vec, file = "labels_vec_for_distmat_20_clusters.RData")
adj_mat_clusters <- mat.new.new[unlist(clusters_groups), unlist(clusters_groups)]
save(adj_mat_clusters, file = 'adj_mat_parsed_clusters_20_clusters.RData')













