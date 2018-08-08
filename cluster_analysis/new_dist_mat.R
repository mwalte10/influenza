#SETUP
{
setwd('~/Desktop/CSHL/data')
load('DNA_GEN.RData')
load('AA_GEN.RData')
library(vwr)
library(igraph)
library(RColorBrewer)
library(DiagrammeR)
library(vegan)
years_vec <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
#load('dist_mat_twenty_clusters.RData')
#load("years.colors.dist_mat.RData")
setwd('~/Desktop/CSHL/cluster_analysis')
load('labels_for_simp.adj.mat.RData')
load("labels_vec_for_distmat_20_clusters.RData")
}

#can all be skipped if everything is ok
{###build up dna_to
{
dna_to <- c()
for(i in 1:length(labels_vec)){
  seq <- unlist(dna_gen[which(dna_gen[,4] == labels_vec[i]),5])[1: (566 * 3)]
  if(any(is.na(seq))){next}
  dna_to[i] <- paste(seq, collapse = '')
  print(i / length(labels_vec))
}
}

##recreating distance_matrix
{
  ham_dist <- function(k){
    comp <- dna_to[k]
    to <- dna_to[k:length(dna_to)]
    ham <- hamming.distance(as.character(comp), as.character(to))
    if(length(ham) == 0){return(0)}
    dist <- unname(ham)
    return(dist)
  }
  
  dist_mat.large <- matrix(NA, ncol = length(dna_to), nrow = length(dna_to))
  for(k in 1:ncol(dist_mat.large)){
    output <- ham_dist(k)
    dist_mat.large[(k:length(dna_to)),k] <- output
    dist_mat.large[k,(k:length(dna_to))] <- output
    print((length(dna_to) - k)^2 / (length(dna_to))^2)
  }

save(dist_mat, file = 'dist_mat_large.RData')

}

#get years.colors
{
years <- c()
for(i in 1:length(labels_vec)){
  years[i] <- dna_gen[which(dna_gen[,4] == labels_vec[i]), 1]
  print(i / length(labels_vec))
}
years <- unlist(years)
years.colors <- c()
colors <- c(brewer.pal(8, "Accent"), "red")
for(i in 1:length(years)){
  years.colors[i] <- colors[which(years_vec == years[i])]
  remove(i)
}
}

##full mds plot
{
counts <- c()
for(i in 1:length(unique(years.colors))){
  counts[i] <- length(which(years.colors == unique(years.colors)[i]))
}
counts.indexes <- list()
for(i in 1:length(unique(years.colors))){
  counts.indexes[[i]] <- which(years.colors == unique(years.colors)[i])
}
weights <- list()
for(i in 1:length(counts.indexes)){
  weights[[i]] <- rep((1/ length(counts.indexes[[i]])), length(counts.indexes[[i]]))
}
weights <- unlist(weights)
weighted_mds <- wcmdscale(dist_mat.large, w = weights, k = 2)
save(weighted_mds, file = "weighted_mds_mat.RData")
load('weighted_mds_mat.RData')

##get cluster formatting
edges <- list()
for(i in 1:ncol(dist_mat.large)){
	right <- which(dist_mat.large[,i] == 1)
	left <- rep(i, length(right))
	edges[[i]] <- cbind(left, right)
}
edges <- do.call(rbind, edges)
# remove <- c(unique(which(weighted_mds[,2] < -200)), unique(which(weighted_mds[,2] > 200)))
# remove_left <- list()
# remove_right <- list()
# for(i in 1:length(remove)){
#   remove_left[[i]] <- which(edges[,1] == remove[i])
#   remove_right[[i]] <- which(edges[,2] == remove[i])
# }
# remove_edges <- unique(c(unlist(remove_left), unlist(remove_right)))
# edges <- edges[-remove_edges,]
vertices <- unique(c(edges[,1], edges[,2]))
graph.object <- graph.data.frame(edges, vertices, directed = FALSE)
layout =  weighted_mds


# layout <- layout[-remove,]
V(graph.object)$name = NA
V(graph.object)$color = years.colors
V(graph.object)$frame.color = 'black'
E(graph.object)$curved = 0
V(graph.object)$size = rep(2, length(vertices))
par(mar=c(0, 0, 2, 0))
plot(graph.object, layout = layout, main = "Top Twenty Clusters \nMultidimensional Scaling", xaxs="i", yaxs="i")
legend(0.6882139, -0.4903915, legend = years_vec, fill = colors, cex = 0.8)




}
  
  
###find the top twenty sequences, highlight with black dots
{{
dates <- as.Date(paste(dna_gen[,1], '-', dna_gen[,2], '-', dna_gen[,3], sep =''), '%Y-%m-%d')
years_twenty <- dna_gen[,1]
which.labels <- c()
for(i in 1:length(labels.new.new)){
  which.labels[i] <- which(dna_gen[,4] == labels.new.new[i])
  print(i / length(labels.new.new))
}
dates <- dates[which.labels]
years_twenty <- years_twenty[which.labels]
years_twenty <- unlist(years_twenty)
ones <- c()
for(i in 1:ncol(dist_mat.large)){
  ones[i] <- length(which(dist_mat.large[,i] == 1))
}
top_twenty_values <- unique(tail(sort(ones), 20))
top_twenty_which <- list()
for(i in 1:length(top_twenty_values)){
  top_twenty_which[[i]] <- which(ones == top_twenty_values[i])
}
top_twenty_which <- unlist(top_twenty_which)
top_twenty_labels <- c()
for(i in 1:length(top_twenty_which)){
  top_twenty_labels[i] <- labels.new.new[top_twenty_which[i]]
}
top_twenty_nuc.main <- c()
top_twenty_aa.main <- c()
for(i in 1:length(top_twenty_which)){
  top_twenty_nuc.main[i] <- dna_gen[which(dna_gen[,4] == top_twenty_labels[i]), 5]
  top_twenty_aa.main[i] <- aa_gen[which(aa_gen[,4] == top_twenty_labels[i]), 5]
}
top_twenty_list <- list()
for(i in 1:length(top_twenty_which)){
  top_twenty_list[[i]] <- labels.new.new[which(dist_mat.large[,top_twenty_which[i]] == 1)]
}
top_twenty_nuc.cluster <- list()
top_twenty_aa.cluster <- list()
for(i in 1:length(top_twenty_which)){
  cluster <- top_twenty_list[[i]]
  vec.nuc <- c()
  vec.aa <- c()
  for(j in 1:length(cluster)){
   vec.nuc[j] <- dna_gen[which(dna_gen[,4] == cluster[j]), 5]
   vec.aa[j] <- aa_gen[which(aa_gen[,4] == cluster[j]), 5]
  }
  top_twenty_nuc.cluster[[i]] <- vec.nuc
  top_twenty_aa.cluster[[i]] <- vec.aa
}
s_ns <- list()
for(i in 1:length(top_twenty_which)){
  vec <- c()
  aa_seqs <- top_twenty_aa.cluster[[i]]
  for(j in 1:length(aa_seqs)){
    vec[j] <- paste(unlist(aa_seqs[j]), collapse = '') == paste(unlist(top_twenty_aa.main[i]), collapse = '')
  }
  s_ns[[i]] <- vec
}
synon_counts <- c()
ns_counts <- c()
for(i in 1:length(top_twenty_which)){
  synon_counts[i] <- length(which(s_ns[[i]] == TRUE))
  ns_counts[i] <- length(which(s_ns[[i]] == FALSE))
}


plot(x = synon_counts, y = ns_counts, pch = 16, main = "Nonsynomous Counts vs. Synonymous Counts \nAll Data Included",
     xlab = "Number of Synonymous Variants", ylab = "Number of Nonsynonymous Variants")
legend(x = 1522.99, y = 461.8832, legend = "y = 159.10556 + 0.08798 * x \nR^2 = 0.3007", cex = 0.8)
x <- seq(1,4000)
lines(x, y = 159.10556 + 0.08798 * x, col = 'red')

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(x = dates[top_twenty_which], y = (ns_counts / synon_counts), pch = 16, xlab = "Date", ylim = c(0, max(ns_counts/synon_counts)),
     ylab = "Nonsynonymous to \nSynonymous Fraction", main = "Nonsynonymous / Synonymous Fraction")
text(dates[top_twenty_which], (ns_counts / synon_counts), labels = labels.new.new[top_twenty_which], pos = 1, cex = 0.5)
legend(x = 96.61695, y = 0.6687728, legend = "Correlation Coefficient: 0.086561", cex = 0.8)

}
  
  top_twenty_aa <- c()
  for(i in 1:length(top_twenty_aa.main)){
    top_twenty_aa[i] <- paste(unlist(top_twenty_aa.main[[i]]), collapse = '')
  }
  top_twenty_aa <- top_twenty_aa[order(dates[top_twenty_which])]
  
  ham_dist_aa <- function(k){
    comp <- top_twenty_aa[k]
    to <- top_twenty_aa[k:length(top_twenty_aa)]
    ham <- hamming.distance(as.character(comp), as.character(to))
    if(length(ham) == 0){return(0)}
    dist <- unname(ham)
    return(dist)
  }
  
  dist_mat_aa <- matrix(NA, ncol = length(top_twenty_aa), nrow = length(top_twenty_aa))
  for(k in 1:ncol(dist_mat_aa)){
    output <- ham_dist_aa(k)
    dist_mat_aa[(k:length(top_twenty_aa)),k] <- output
    dist_mat_aa[k,(k:length(top_twenty_aa))] <- output
    print((length(top_twenty_aa) - k)^2 / (length(top_twenty_aa))^2)
  }
  # ns_counts <- ns_counts[order(dates[top_twenty_which])]
  # synon_counts <- synon_counts[order(dates[top_twenty_which])]
  fractions <- ns_counts / synon_counts
  fractions <- fractions[order(dates[top_twenty_which])]
  labels.twenty <- labels.new.new[top_twenty_which]
  labels.twenty <- labels.twenty[order(dates[top_twenty_which])]
   synon_counts <- synon_counts[order(dates[top_twenty_which])]
   ns_counts <- ns_counts[order(dates[top_twenty_which])]
   dates <- dates[top_twenty_which]
   dates <- dates[order(top_twenty_which)]
   years_twenty <- years_twenty[top_twenty_which]
   years_twenty <- years_twenty[order(top_twenty_which)]
   years_twenty.colors <- c()
   for(i in 1:length(years_twenty)){
    years_twenty.colors[i] <- colors[which(years_twenty[i] == years_vec)]
   }
   
{  par(mfrow = c(3,1), xpd = FALSE)
  plot(c(fractions[5], fractions[6], fractions[7], fractions[8]), 
       xaxt = 'n', xlab = NA, col = c( years_twenty.colors[5],  years_twenty.colors[6],
                                       years_twenty.colors[7],  years_twenty.colors[8])
       ylim = c(0,1), ylab = 'NS / S', pch = 16)
  text(seq(1,4), c(fractions[5], fractions[6], fractions[7], fractions[8]), 
       labels = c(labels.twenty[5], labels.twenty[6], labels.twenty[7], labels.twenty[8]), 
       pos = 3, cex = 0.6)
  text(seq(1,3), c(fractions[5], fractions[6], fractions[7]),
       labels = c('p = 0.9469', 'p = 0.3104', 'p = 0.07601'), 
       pos = 1, cex = 0.6)
  abline(v = 3.5)
  plot(c(fractions[10], fractions[16], fractions[11]), 
       xaxt = 'n', xlab = NA, 
       ylim = c(0,1), ylab = 'NS / S', pch = 16)
  text(seq(1,3), c(fractions[10], fractions[16], fractions[11]), 
       labels = c(labels.twenty[10], labels.twenty[16], labels.twenty[11]), 
       pos = 3, cex = 0.6)
  text(seq(1,2),  c(fractions[10], fractions[16]),
       labels = c('p = 0.04358', 'p = 0.001203'), 
       pos = 1, cex = 0.6)
  abline(v = 2.5)
  plot(c(fractions[13], fractions[21], fractions[18], fractions[22]), 
       xaxt = 'n', xlab = NA,
       ylim = c(0,1), ylab = 'NS / S', pch = 16)
  text(seq(1,3), c(fractions[13], fractions[21], fractions[18]), 
       labels = c('p = 0.1007', 'p = 0.4926', 'p = 0.0761'), 
       pos = 1, cex = 0.6)
  text(seq(1,4), c(fractions[13], fractions[21], fractions[18], fractions[22]), 
       labels = c(labels.twenty[13], labels.twenty[21], labels.twenty[18], labels.twenty[22]), 
       pos = 3, cex = 0.6)
  abline(v = 3.5)}
   par(mar = mar.default + c(0, 0, 4, 0))
  par(mfrow = c(2,2), xpd = TRUE, mar = mar.default + c(0, 4, 0, 0))
 ##plot A 
{  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0)) 
  plot(x = c(1, 1, 1, 2, 5, 5, 5, 6),
    y = c(fractions[5], fractions[6], fractions[7], fractions[8], fractions[13], fractions[18], fractions[21], fractions[22]), 
       xaxt = 'n', xlab = 'Hamming Distance', 
    col = c(years_twenty.colors[5], years_twenty.colors[6], years_twenty.colors[7],
            years_twenty.colors[8], years_twenty.colors[13], years_twenty.colors[18],
            years_twenty.colors[21], years_twenty.colors[22]),
       ylim = c(0,1), ylab = 'Nonsynonymous to synonymous ratio \n(NS / S)', pch = 16)
  ##add conf int
  ci_5 <- prop.test(sum(ns_counts[5], ns_counts[6], ns_counts[7]),sum(synon_counts[5], synon_counts[7], synon_counts[9]), correct = FALSE)
  arrows(1, ci_5$conf.int[1], 1, ci_5$conf.int[2], length=0.05, angle=90, code=3)
  ci_8 <- prop.test(c(ns_counts[8]),c(synon_counts[8]), correct = FALSE)
  arrows(2, ci_8$conf.int[1], 2, ci_8$conf.int[2], length=0.05, angle=90, code=3)
  ci_13 <- prop.test(sum(ns_counts[13], ns_counts[18], ns_counts[21]),sum(synon_counts[13], synon_counts[18], synon_counts[21]), correct = FALSE)
  arrows(5, ci_13$conf.int[1], 5, ci_13$conf.int[2], length=0.05, angle=90, code=3)
  ci_22 <- prop.test(c(ns_counts[22]),c(synon_counts[22]), correct = FALSE)
  arrows(6, ci_22$conf.int[1], 6, ci_22$conf.int[2], length=0.05, angle=90, code=3)
  
  ##add in distance bars
  arrows(1, 0 , 1.98, 0, angle = 90, length=0.05, code = 3)
  arrows(2, 0 , 4.98, 0, angle = 90, length=0.05, code = 3)
  arrows(5, 0 , 5.98, 0, angle = 90, length=0.05, code = 3)
  text(x= 1.5, y = 0.03, 'Distance = 1', cex = 0.6)
  text(x = 3.5, y = 0.08, '*')
  text(x= 3.5, y = 0.03, 'Distance = 3', cex = 0.6)
  text(x= 5.5, y = 0.03, 'Distance = 1', cex = 0.6)
  
  
  text(x = c(1, 1, 1, 2), y = c(fractions[5], fractions[6], fractions[7], fractions[8]), 
       labels = c(labels.twenty[5], labels.twenty[6], labels.twenty[7], labels.twenty[8]), 
       col = c(years_twenty.colors[5], years_twenty.colors[6], years_twenty.colors[7],
               years_twenty.colors[8]),
       pos = 4, cex = 0.5, font = 2)
  text(x= c(5,5,5,6), y = c( fractions[13], fractions[18], fractions[21], fractions[22]), c(labels.twenty[13], labels.twenty[18], labels.twenty[21], labels.twenty[22]), 
       col = c(years_twenty.colors[13],years_twenty.colors[28], years_twenty.colors[21], years_twenty.colors[22]), 
       pos = 2, cex = 0.5, font = 2)
  text(x = c(1.5, 3.5, 5.5), y = rep(-0.13, 3), 
       labels = c("r261l", "t131k \nr142k \nl261q", "a212t"), 
       cex = 0.6, col = 'black')}

 ##plot B 
 { mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0)) 
  plot(x = c(1, 4, 7),
       c(fractions[9], fractions[15], fractions[17]), 
       xaxt = 'n', xlab = 'Hamming Distance', 
       col = c(years_twenty.colors[9], years_twenty.colors[15], years_twenty.colors[17]),
       ylim = c(0,1), ylab = 'Nonsynonymous to synonymous ratio \n(NS / S)', pch = 16)
  ##add conf int
  ci_9 <- prop.test(c(ns_counts[9]),c(synon_counts[9]), correct = FALSE)
  arrows(1, ci_9$conf.int[1], 1, ci_9$conf.int[2], length=0.05, angle=90, code=3)
  ci_15 <- prop.test(c(ns_counts[15]),c(synon_counts[15]), correct = FALSE)
  arrows(4, ci_15$conf.int[1], 4, ci_15$conf.int[2], length=0.05, angle=90, code=3)
  ci_17 <- prop.test(c(ns_counts[17]),c(synon_counts[17]), correct = FALSE)
  arrows(7, ci_17$conf.int[1], 7, ci_17$conf.int[2], length=0.05, angle=90, code=3)
  
  arrows(1, 0 , 3.98, 0, angle = 90, length=0.05, code = 3)
  arrows(4, 0 , 6.98, 0, angle = 90, length=0.05, code = 3)
  text(x= 2.5, y = 0.03, 'Distance = 3', cex = 0.6)
  text(x = 2.5, y = 0.08, '***')
  text(x= 5.5, y = 0.03, 'Distance = 3', cex = 0.6)

  text(x = c(1, 4), c(fractions[9], fractions[15]), 
       labels = c(labels.twenty[9], labels.twenty[15]), 
       col = c(years_twenty.colors[9], years_twenty.colors[15]),
       pos = 4, cex = 0.5, font = 2)
  text(7, fractions[17], labels = labels.twenty[17], col = years_twenty.colors[17], pos = 2, cex = 0.5, font = 2)
  text(x = c(2.5, 5.5), y = rep(-0.13, 2), 
       labels = c("n121k \ni140m \ng479e", "s46t \nt135k \nm140i"), 
       cex = 0.6, col = 'black')}

#plot c  
{  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0)) 
  plot(x = c(1, 4, 7),
       c(fractions[9], fractions[12], fractions[19]), 
       col = c(years_twenty.colors[9], years_twenty.colors[12], years_twenty.colors[19]),
       xaxt = 'n', xlab = 'Hamming Distance', 
       ylim = c(0,1), ylab = 'Nonsynonymous to synonymous ratio \n(NS / S)', pch = 16)
  ##add conf int
  arrows(1, ci_9$conf.int[1], 1, ci_9$conf.int[2], length=0.05, angle=90, code=3)
  ci_12 <- prop.test(c(ns_counts[12]),c(synon_counts[12]), correct = FALSE)
  arrows(4, ci_12$conf.int[1], 4, ci_12$conf.int[2], length=0.05, angle=90, code=3)
  ci_19 <- prop.test(c(ns_counts[19]),c(synon_counts[19]), correct = FALSE)
  arrows(7, ci_19$conf.int[1], 7, ci_19$conf.int[2], length=0.05, angle=90, code=3)
  
  arrows(1, 0 , 3.98, 0, angle = 90, length=0.05, code = 3)
  arrows(4, 0 , 6.98, 0, angle = 90, length=0.05, code = 3)
  text(x= 2.5, y = 0.03, 'Distance = 3', cex = 0.6)
  text(x = 2.5, y = 0.08, '**')
  text(x= 5.5, y = 0.03, 'Distance = 3', cex = 0.6)
  text(x = 5.5, y = 0.08, '***')
  
  text(x = c(1, 4), c(fractions[9], fractions[12]), 
       labels = c(labels.twenty[9], labels.twenty[12]), 
       col = c(years_twenty.colors[9], years_twenty.colors[12]), font = 2,
       pos = 4, cex = 0.5)
  text(7, fractions[19], labels = labels.twenty[19], col = years_twenty.colors[19],
       font = 2, pos = 2, cex = 0.5)
  text(x = c(2.5, 5.5), y = rep(-0.13, 2), 
       labels = c("k92r \nn121k \nh311q", "e62g \nt135k \nr142g"), 
       cex = 0.6, col = 'black')}

#plot d
{  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 2, 0, 0)) 
  plot(x = c(1, 2, 3, 3),
       y = c(fractions[9], fractions[11], fractions[10], fractions[16]), 
       xaxt = 'n', xlab = 'Hamming Distance',
       col = c(years_twenty.colors[9], years_twenty.colors[11], years_twenty.colors[10], years_twenty.colors[16]),
       ylim = c(0,1), ylab = 'Nonsynonymous to synonymous ratio \n(NS:S)', pch = 16)
  ##add conf int
  arrows(1, ci_9$conf.int[1], 1, ci_9$conf.int[2], length=0.05, angle=90, code=3)
  ci_10 <- prop.test(sum(ns_counts[10], ns_counts[16]),sum(synon_counts[10], synon_counts[16]), correct = FALSE)
  arrows(3, ci_10$conf.int[1], 3, ci_10$conf.int[2], length=0.05, angle=90, code=3)
  ci_11 <- prop.test(c(ns_counts[11]),c(synon_counts[11]), correct = FALSE)
  arrows(2, ci_11$conf.int[1], 2, ci_11$conf.int[2], length=0.05, angle=90, code=3)
  
  arrows(1, 0 , 1.98, 0, angle = 90, length=0.05, code = 3)
  arrows(2, 0 , 2.98, 0, angle = 90, length=0.05, code = 3)
  text(x= 1.5, y = 0.03, 'Distance = 1', cex = 0.6)
  text(x = 1.5, y = 0.08, '***')
  text(x= 2.5, y = 0.03, 'Distance = 1', cex = 0.6)
  text(x = 2.5, y = 0.08, '**')
  
  text(x = c(1, 2), c(fractions[9], fractions[11]), 
       labels = c(labels.twenty[9],  labels.twenty[11]), 
       col = c(years_twenty.colors[9], years_twenty.colors[11]),
       pos = 4, cex = 0.5, font =2)
  text(c(3,3), c(fractions[10], fractions[16]), labels = c(labels.twenty[10], labels.twenty[16]),
       col = c(years_twenty.colors[10], years_twenty.colors[16]),
       pos = 2, cex = 0.5, font =2)
  text(x = c(1.5, 2.5), y = rep(-0.13, 2), 
       labels = c("r142g", "n121k"), 
       cex = 0.6, col = 'black')
  
  }
  legend(-.2,-0.3,inset = 0,
         legend = years_vec, 
         col=colors, lwd=5, cex=.4, horiz = TRUE,xpd="NA",
         title = "Years")
  
  #two sample binomial tests
  {
    ##first graph
    prop.test(c(sum(ns_counts[5], ns_counts[6], ns_counts[7]),ns_counts[8]),
              c(sum(synon_counts[5], synon_counts[6],synon_counts[7]),synon_counts[8]),correct=FALSE)

    prop.test(c(sum(ns_counts[13], ns_counts[18], ns_counts[21]),ns_counts[8]),
              c(sum(synon_counts[13], synon_counts[18],synon_counts[21]),synon_counts[8]),correct=FALSE)
    prop.test(c(sum(ns_counts[13], ns_counts[18], ns_counts[21]),ns_counts[22]),
              c(sum(synon_counts[13], synon_counts[18],synon_counts[21]),synon_counts[22]),correct=FALSE)

    
    prop.test(c(ns_counts[9],ns_counts[15]),c(synon_counts[9],synon_counts[15]),correct=FALSE)
    prop.test(c(ns_counts[15],ns_counts[17]),c(synon_counts[15],synon_counts[17]),correct=FALSE)
    ##example of CI
    ci_9 <- add4ci(c(ns_counts[9]),c(synon_counts[9]), 0.95)
    arrows(fractions[9], ci_9$conf.int[1], fractions[9],ci_9$conf.int[2], length=0.05, angle=90, code=3)
    
    prop.test(c(ns_counts[9],ns_counts[12]),c(synon_counts[9],synon_counts[12]),correct=FALSE)
    prop.test(c(ns_counts[12],ns_counts[19]),c(synon_counts[12],synon_counts[19]),correct=FALSE)
    
    prop.test(c(ns_counts[11],sum(ns_counts[10], ns_counts[16])),c(synon_counts[11],sum(synon_counts[10], synon_counts[16])),correct=FALSE)
    prop.test(c(ns_counts[11],ns_counts[9]),c(synon_counts[11],synon_counts[9]),correct=FALSE)
    
    
  }
  
  
  for(i in 1:566){
    from <- unlist(strsplit(top_twenty_aa[12], split = ''))
    to <- unlist(strsplit(top_twenty_aa[19], split = ''))
    if(from[i] == to[i]){next}else{print(c(from[i], (i), to[i]))}
  }
  }
  
#Look for K at position 142
{
aa_seq_list <- list()
for(i in 1:length(labels_vec)){
  aa_seq_list[i] <- aa_gen[which(aa_gen[,4] == labels_vec[i]), 5]
}
seq_11.edges <- which(edges[,1] == 5233)
edges_11.move <- edges[seq_11.edges,]
edges <- edges[-seq_11.edges,]
vertices <- vertices[-5233]
layout_11.move <- layout[5233,]
layout <- layout[-5233,]
aa_seq_list <- aa_seq_list[-5233]

# pos.142 <- c()
# for(i in 1:length(aa_seq_list)){
#   pos.142[i] <- aa_seq_list[[i]][142 + 16]
# }
# list.142 <- list()
# for(i in 1:length(which(pos.142 == 'g'))){
#   list.142[[i]] <- which(edges[,1] == which(pos.142 == 'g')[i])
# }
# list.142 <- unlist(list.142)
pos.137 <- c()
for(i in 1:length(aa_seq_list)){
  pos.137[i] <- aa_seq_list[[i]][121 + 16]
}
list <- list()
for(i in 1:length(which(pos.137 == 'k'))){
  list[[i]] <- which(edges[,1] == which(pos.137 == 'k')[i])
}
list <- unlist(list)

edges.move <- edges[list,]
edges <- edges[-list,]
edges <- rbind(edges, edges.move)
edges <- rbind(edges, edges_11.move)
#vertices <- unique(c(edges[,1], edges[,2]))
vertices.move <- vertices[which(pos.137 == 'k')]
vertices <- vertices[-which(pos.137 == 'k')]
vertices <- c(vertices, vertices.move)
vertices <- c(vertices, 5233)
graph.object <- graph.data.frame(edges, vertices, directed = FALSE)
#layout =  weighted_mds[,c(1,2)]
layout.move <- layout[which(pos.137 == 'k'),]
layout <- layout[-which(pos.137 == 'k'),]
layout <- rbind(layout, layout.move)
layout <- rbind(layout, layout_11.move)
V(graph.object)$name = NA
V(graph.object)$color = c(rep('black', (length(vertices) - length(vertices.move) - 1)), rep('deepskyblue3', length(vertices.move)), "darkseagreen1")
V(graph.object)$frame.color = 'black'
E(graph.object)$curved = 0
V(graph.object)$size = c(rep(3, (length(vertices) - 1)), 5)
plot(graph.object, layout = layout, main = "Top Twenty Clusters \nk at Position 121")
text(x= 0.1, y = -0.3673939, "EPI_ISL_225827", font = 2)
legend(0.6152743, -0.9535271 , legend = c('k at position 121'), fill = c('deepskyblue3') , cex = 0.8)

       }

#highlight clade a1b and a2
  which_a1a <- c()
  which_a1a_edges <- list()
  for(i in 1:length(clade_present_a1a)){
    which_a1a[i] <- which(labels_vec == clade_present_a1a[i]) 
    which_a1a_edges[[i]] <- which(edges[,1] == which_a1a[i])
  }
  which_a1a_edges <- unlist(which_a1a_edges)
  
which_a1b <- c()
which_a1b_edges <- list()
for(i in 1:length(clade_present_a1b)){
  which_a1b[i] <- which(labels_vec == clade_present_a1b[i]) 
  which_a1b_edges[[i]] <- which(edges[,1] == which_a1b[i])
}
which_a1b_edges <- unlist(which_a1b_edges)

which_a2 <- c()
which_a2_edges <- list()
for(i in 1:length(clade_present_a2)){
  which_a2[i] <- which(labels_vec == clade_present_a2[i]) 
  which_a2_edges[[i]] <- which(edges[,1] == which_a2[i])
}
which_a2_edges <- unlist(which_a2_edges)

which_a3 <- c()
which_a3_edges <- list()
for(i in 1:length(clade_present_a3)){
  which_a3[i] <- which(labels_vec == clade_present_a3[i]) 
  which_a3_edges[[i]] <- which(edges[,1] == which_a3[i])
}
which_a3_edges <- unlist(which_a3_edges)

which_a4 <- c()
which_a4_edges <- list()
for(i in 1:length(clade_present_a3)){
  which_a4[i] <- which(labels_vec == clade_present_a4[i]) 
  which_a4_edges[[i]] <- which(edges[,1] == which_a4[i])
}
which_a4_edges <- unlist(which_a4_edges)

# move_edges <- c(which_a1a_edges, which_a1b_edges, which_a2_edges, which_a3_edges, which_a4_edges)
# move_vertices <- unique(c(which_a1a, which_a1b, which_a2, which_a3, which_a4))
# move_edges <- c()
# for(i in 1:length(move_vertices)){
#   move_edges[i] <- which(edges[,1] == move_vertices[i])
# }

move_edges <- c(which_a1b_edges, which_a2_edges)
move_vertices <- unique(c(which_a1b, which_a2))
move_edges <- c()
for(i in 1:length(move_vertices)){
  move_edges[i] <- which(edges[,1] == move_vertices[i])
}


colors_clusters <- c(brewer.pal(5, "Accent"))
edges.move <- edges[move_edges,]
edges <- edges[-move_edges,]
edges <- rbind(edges, edges.move)
vertices.move <- vertices[move_vertices]
vertices <- vertices[-move_vertices]
vertices <- c(vertices, vertices.move)
graph.object <- graph.data.frame(edges, vertices, directed = FALSE)
layout.move = layout[move_vertices,]
layout = layout[-move_vertices,]
layout = rbind(layout, layout.move)
V(graph.object)$name = NA
# V(graph.object)$color = c(rep('black', (length(vertices) - length(move_vertices))), rep(colors_clusters[1], length(which_a1a)),
#                           rep(colors_clusters[2], length(which_a1b)), rep(colors_clusters[3], length(which_a2)),
#                           rep(colors_clusters[4], length(which_a3)), rep(colors_clusters[5], length(which_a4)))
V(graph.object)$color = c(rep('black', (length(vertices) - length(move_vertices))), 
                          rep('seagreen', length(which_a1b)), rep('deepskyblue3', length(which_a2)))
V(graph.object)$frame.color = 'black'
E(graph.object)$curved = 0
V(graph.object)$size = rep(3, length(vertices))
plot(graph.object, layout = layout, main = "Top Twenty Clusters \nMultidimensional Scaling")
legend('bottomright', legend = c('clade 3ca.a1b', 'clade 3ca.a2'), fill = c('seagreen', 'deepskyblue3'))

  
  
{  ##Subclade 3ca.a1a
  {site_135_a1a <- c()
    for(i in 1:nrow(aa_gen)){
      site_135_a1a[i] <- unlist(aa_gen[i,5])[135 + 16]
    }
    labels_clade_a1a <- unlist(aa_gen[which(site_135_a1a == 'k'),4])
  }
  
  ##Subclade 3ca.a1b
  {site_92_a1b <- c()
    for(i in 1:nrow(aa_gen)){
      site_92_a1b[i] <- unlist(aa_gen[i,5])[92 + 16]
    }
    labels_92_a1b <- unlist(aa_gen[which(site_92_a1b == 'r'),4])
    
    site_311_a1b <- c()
    for(i in 1:nrow(aa_gen)){
      site_311_a1b[i] <- unlist(aa_gen[i,5])[311 + 16]
    }
    labels_311_a1b <- unlist(aa_gen[which(site_311_a1b == 'q'),4])
    
    labels_clade_a1b <- Reduce(intersect, list(labels_311_a1b,labels_92_a1b))
    
  }
  
  ##Subclade 3c2.a2
  {site_131_a2 <- c()
    for(i in 1:nrow(aa_gen)){
      site_131_a2[i] <- unlist(aa_gen[i,5])[131 + 16]
    }
    labels_131_a2 <- unlist(aa_gen[which(site_131_a2 == 'k'),4])
    
    
    site_142_a2 <- c()
    for(i in 1:nrow(aa_gen)){
      site_142_a2[i] <- unlist(aa_gen[i,5])[142 + 16]
    }
    labels_142_a2 <- unlist(aa_gen[which(site_142_a2 == 'k'),4])
    
    site_261_a2 <- c()
    for(i in 1:nrow(aa_gen)){
      site_261_a2[i] <- unlist(aa_gen[i,5])[261 + 16]
    }
    labels_261_a2 <- unlist(aa_gen[which(site_261_a2 == 'q'),4])
    
    labels_clade_a2 <- Reduce(intersect, list(labels_131_a2,labels_142_a2,labels_261_a2))
    
  }
  
  ##Subclade 3c2.a3
  {site_121_a3 <- c()
    for(i in 1:nrow(aa_gen)){
      site_121_a3[i] <- unlist(aa_gen[i,5])[121 + 16]
    }
    labels_121_a3 <- unlist(aa_gen[which(site_121_a3 == 'k'),4])
    
    site_144_a3 <- c()
    for(i in 1:nrow(aa_gen)){
      site_144_a3[i] <- unlist(aa_gen[i,5])[144 + 16]
    }
    labels_144_a3 <- unlist(aa_gen[which(site_144_a3 == 'k'),4])
    
    site_135_a3 <- c()
    for(i in 1:nrow(aa_gen)){
      site_135_a3[i] <- unlist(aa_gen[i,5])[135 + 16]
    }
    labels_135_a3 <- unlist(aa_gen[which(site_135_a3 == 'k'),4])
    
    labels_clade_a3 <-  Reduce(intersect, list(labels_121_a3,labels_144_a3,labels_135_a3))
    
  }
  
  ##Subclade 3c2.a4
  {site_53_a4 <- c()
    for(i in 1:nrow(aa_gen)){
      site_53_a4[i] <- unlist(aa_gen[i,5])[53 + 16]
    }
    labels_53_a4 <- unlist(aa_gen[which(site_53_a4 == 'n'),4])
    
    site_144_a4 <- c()
    for(i in 1:nrow(aa_gen)){
      site_144_a4[i] <- unlist(aa_gen[i,5])[144 + 16]
    }
    labels_144_a4 <- unlist(aa_gen[which(site_144_a4 == 'r'),4])
    
    site_171_a4 <- c()
    for(i in 1:nrow(aa_gen)){
      site_171_a4[i] <- unlist(aa_gen[i,5])[171 + 16]
    }
    labels_171_a4 <- unlist(aa_gen[which(site_171_a4 == 'k'),4])
    
    site_192_a4 <- c()
    for(i in 1:nrow(aa_gen)){
      site_192_a4[i] <- unlist(aa_gen[i,5])[192 + 16]
    }
    labels_192_a4 <- unlist(aa_gen[which(site_192_a4 == 't'),4])
    
    site_197_a4 <- c()
    for(i in 1:nrow(aa_gen)){
      site_197_a4[i] <- unlist(aa_gen[i,5])[197 + 16]
    }
    labels_197_a4 <- unlist(aa_gen[which(site_197_a4 == 'h'),4])
    
    labels_clade_a4 <- Reduce(intersect, list(labels_53_a4, labels_144_a4, labels_171_a4, labels_197_a4))
    
  }
  
  clade_3c2.a <- list('3ca.a1a' = labels_clade_a1a, '3ca.a1b' = labels_clade_a1b, '3c2.a2' = labels_clade_a2, 
                      '3c2.a3' = labels_clade_a3, '3c2.a4' = labels_clade_a4)
  
  clade_present_a1a <- Reduce(intersect, list(labels_vec[vertices], clade_3c2.a$'3ca.a1a'))
  clade_present_a1b <- Reduce(intersect, list(labels_vec[vertices], clade_3c2.a$'3ca.a1b'))
  clade_present_a2 <- Reduce(intersect, list(labels_vec[vertices], clade_3c2.a$'3c2.a2'))
  clade_present_a3 <- Reduce(intersect, list(labels_vec[vertices], clade_3c2.a$'3c2.a3'))
  clade_present_a4 <- Reduce(intersect, list(labels_vec[vertices], clade_3c2.a$'3c2.a4'))
 } 
  
  
{###get all amino acid sequences
clustering <- unname(components(graph.object)$membership)
for(i in 1:20){
  assign(paste('labels_cluster_', i, sep = ''), labels_vec[which(clustering == i)])
  assign(paste('sequence_cluster_', i, sep = ''), labels_vec[which(clustering == i)])
}


##testing to find rates of NS/ S mutation within the clusters
{{
{f <- c('ttt', 'ttc', rep(NA, 4))
l <- c('tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg')
s <- c('tct', 'tcc', 'tca', 'tcg', 'agt', 'agc')
y <- c('tat', 'tac', rep(NA, 4))
c <- c('tgt', 'tgc', rep(NA, 4))
w <- c('tgg', rep(NA, 5))
p <- c('cct', 'ccc', 'cca', 'ccg', rep(NA, 2))
h <- c('cat', 'cac', rep(NA, 4))
q <- c('caa', 'cag', rep(NA, 4))
r <- c('cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg')
i <- c('att', 'atc', 'ata', rep(NA, 3))
m <- c('atg', rep(NA, 5))
t <- c('act', 'acc', 'aca', 'acg', rep(NA, 2))
n <- c('aat', 'aac', rep(NA, 4))
k <- c('aaa', 'aag', rep(NA, 4))
v <- c('gtt', 'gtc', 'gta', 'gtg', rep(NA, 2))
a <- c('gct', 'gcc', 'gca', 'gcg', rep(NA, 2))
d <- c('gat', 'gac', rep(NA, 4))
e <- c('gaa', 'gag', rep(NA, 4))
g <- c('ggt', 'ggc', 'gga', 'ggg', rep(NA, 2))
stop <- c('tga', 'taa', 'tag', rep(NA, 3))
amino_acids <- cbind(f, l, s, y,
                     c, w, p, h,
                     q, r, i, m,
                     t, n, k, v,
                     a, d, e, g,
                     stop)
names_aa <- colnames(amino_acids)
names_aa <- names_aa[1:20]
}
  synon <- c(1, 5, 5, 1, 1, 0, 3, 1, 1, 5, 2, 0, 3,
             1, 1, 3, 3, 1, 1, 3)
  names(synon) <- names_aa
  non_synon <- rep(9, 20) - synon
  
  ####Def am not doing this correctly, look back at it
  seq <- c()
  for(i in 1:length(top.labels)){
    seq[i] <- (aa_gen[which(unlist(aa_gen[,4]) == top.labels[i]), 5])
  }

  synon_rate <- rep(0, length(seq))
  for(i in 1:length(seq)){
    seq_this <- unlist(seq[i])
    for(j in 1:length(seq_this)){
      if(seq_this[j] == 'x'){next}
      synon_rate[i] <- (unname(synon[which(names_aa == seq_this[j])]) + synon_rate[i]) 
    }
  }
  synon_rate <- synon_rate / (566 * 9)
  non_synon_rate <- 1 - synon_rate
}
  
all_list <- list()
for(i in 1:20){
  cluster <- list.neighbors[[i]]
  cluster_list <- list()
  for(j in 1:length(cluster)){
    index <- which(labels_vec == cluster_one[j])
    neighbors <- labels_vec[which(dist_mat[,index] == 1)]
    main <- which(aa_gen[,4] == cluster_one[j])
    main_seq <- paste(unlist(aa_gen[main,5]), collapse = '')
    synon <- c()
    for(k in 1:length(neighbors)){
      neighbor <- which(aa_gen[,4] == neighbors[k])
      neighbor_seq <- paste(unlist(aa_gen[neighbor,5]), collapse = '')
      synon[k] <- main_seq == neighbor_seq
    }
    cluster_list[[j]] <- unlist(synon)
  }
  all_list[[i]] <- unlist(cluster_list)
}

obs_synon <- c()
for(i in 1:20){
  obs_synon[i] <- length(which(all_list[[i]] == TRUE)) / length(all_list[[i]])
}

t.test(obs_synon, synon_rate)
}

##test correlation between cluster size and obs_synonymous rate
sizes <- c()
for(i in 1:20){
  sizes[i] <- length(list.neighbors[[i]])
}
plot(x = sizes, y = obs_synon, type = 'l', xlab = "Cluster \nSize", ylab = "Observed Synonymous \nMutation Rate",
     main = "Correlation between cluster size \nand observed mutation rate",
     ylim = c(0.5, 0.9))
points(x = sizes, y = obs_synon, pch = 16)
cor(obs_synon, sizes)

##find years of top labels
top.years <- c()
for(i in 1:length(top.labels)){
  top.years[i] <- dna_gen[which(unlist(dna_gen[,4]) == top.labels[i]), 1]
}
top.years <- unlist(top.years)

##compare NS/S for clusters
top.sequences <- c()
for(i in 1:length(top.labels)){
  top.sequences[i] <- paste(unlist(dna_gen[which(top.labels[i] == dna_gen[,4]),5]), collapse = '')
}
ns_s_ratio <- (1 - obs_synon) / (obs_synon)
dist_mat_top <- dist_mat[top.which, top.which]

##tests whether time is on the x axis
{
counts <- c()
for(i in 1:length(unique(years.colors.new))){
  counts[i] <- length(which(years.colors.new == unique(years.colors.new)[i]))
}
keep <- list()
for(i in 1:length(unique(years.colors.new))){
  keep[[i]] <- sample(which(years.colors.new == unique(years.colors.new)[i]), 80)
}
keep <- unlist(keep)
testing_temp <- new.test[keep, keep]
years.colors.new <- years.colors.new[keep]
mds_temp <- cmdscale(testing_temp)
plot(mds_temp, col = years.colors.new, pch = 16)}


}











