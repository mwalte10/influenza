
{
  nine_neighbors <- labels.new.new[which(dist_mat.large[,which(labels.new.new == labels.twenty[9])] == 1)]
nine_neighbors_seq <- c()
for(i in 1:length(nine_neighbors)){
 nine_neighbors_seq[i] <- aa_gen[which(aa_gen[,4] == nine_neighbors[i]),5]
}
pos.121 <- c()
pos.142 <- c()
for(i in 1:length(nine_neighbors)){
  pos.121[i] <- unlist(nine_neighbors_seq[i])[121 + 16]
  pos.142[i] <- unlist(nine_neighbors_seq[i])[142 + 16]
}

x <- dist_mat.large[,which(nine_neighbors[which(pos.142 == 'g')] == labels.new.new)]
main_connector_label <- labels.new.new[which(nine_neighbors[which(pos.142 == 'g')] == labels.new.new)]
main_seq <- paste(unlist(aa_gen[which(aa_gen[,4] == main_connector_label),5]), collapse = '')
connector_labels <- labels.new.new[which(x == 1)]
connector_seqs <- c()
for(i in 1:length(connector_labels)){
  connector_seqs[i] <- paste(unlist(aa_gen[which(aa_gen[,4] == connector_labels[i]),5]), collapse = '')
}
binary_nss <- c()
for(i in 1:length(connector_labels)){
  binary_nss[i] <- main_seq == connector_seqs[i]
}

connector_two_seq <- c()
for(i in 1:length(binary_nss)){
  connector_two_seq[i] <- aa_gen[which(aa_gen[,4] == connector_labels[i]),5]
}
pos.121 <- c()
for(i in 1:length(connector_two_seq)){
  pos.121[i] <- unlist(connector_two_seq[[i]])[121 + 16]
}
}

###Checking nines neighbors
{
  nine_neighbors <- labels.new.new[which(dist_mat.large[,which(labels.new.new == labels.twenty[9])] == 1)]  
  nine_neighbors_seq <- c()
  for(i in 1:length(nine_neighbors)){
    nine_neighbors_seq[i] <- aa_gen[which(aa_gen[,4] == nine_neighbors[i]),5]
  }
  pos.92 <- c()
  pos.311 <- c()
  pos.121 <- c()
  pos.140 <- c()
  pos.479 <- c()
  for(i in 1:length(nine_neighbors)){
    pos.92[i] <- unlist(nine_neighbors_seq[i])[92 + 16]
    pos.121[i] <- unlist(nine_neighbors_seq[i])[121 + 16]
    pos.140[i] <- unlist(nine_neighbors_seq[i])[140 + 16]
    pos.121[i] <- unlist(nine_neighbors_seq[i])[121 + 16]
    pos.311[i] <- unlist(nine_neighbors_seq[i])[311 + 16]
    pos.479[i] <- unlist(nine_neighbors_seq[i])[479 + 16]
    
  }
  #####STOPPED HERE
  x <- dist_mat.large[,which(nine_neighbors[which(pos.479 == 'e')] == labels.new.new)]
  main_connector_label <- labels.new.new[which(nine_neighbors[which(pos.479 == 'e')] == labels.new.new)]
  main_seq <- paste(unlist(aa_gen[which(aa_gen[,4] == main_connector_label),5]), collapse = '')
  connector_labels <- labels.new.new[which(x == 1)]
  connector_seqs <- c()
  for(i in 1:length(connector_labels)){
    connector_seqs[i] <- paste(unlist(aa_gen[which(aa_gen[,4] == connector_labels[i]),5]), collapse = '')
  }
  binary_nss <- c()
  for(i in 1:length(connector_labels)){
    binary_nss[i] <- main_seq == connector_seqs[i]
  }
  
  connector_two_seq <- c()
  for(i in 1:length(binary_nss)){
    connector_two_seq[i] <- aa_gen[which(aa_gen[,4] == connector_labels[i]),5]
  }
  pos.121 <- c()
  for(i in 1:length(connector_two_seq)){
    pos.121[i] <- unlist(connector_two_seq[[i]])[121 + 16]
  }
}

##checking 12's neighbors
{
  
  twelve_neighbors <- labels.new.new[which(dist_mat.large[,which(labels.new.new == labels.twenty[12])] == 2)]  
  twelve_neighbors_seq <- c()
  for(i in 1:length(twelve_neighbors)){
    twelve_neighbors_seq[i] <- aa_gen[which(aa_gen[,4] == twelve_neighbors[i]),5]
  }
  pos.62 <- c()
  pos.135 <- c()
  pos.142 <- c()
  for(i in 1:length(twelve_neighbors)){
    pos.62[i] <- unlist(twelve_neighbors_seq[i])[62 + 16]
    pos.135[i] <- unlist(twelve_neighbors_seq[i])[135 + 16]
    pos.142[i] <- unlist(twelve_neighbors_seq[i])[142 + 16]
    
  }
  #####STOPPED HERE
  x <- dist_mat.large[,which(twelve_neighbors[which(pos.142 == 'g')[1]] == labels.new.new)]
  main_connector_label <- labels.new.new[which(twelve_neighbors[which(pos.142 == 'g')[1]] == labels.new.new)]
  main_seq <- paste(unlist(aa_gen[which(aa_gen[,4] == main_connector_label),5]), collapse = '')
  connector_labels <- labels.new.new[which(x == 1)]
  connector_seqs <- c()
  for(i in 1:length(connector_labels)){
    connector_seqs[i] <- paste(unlist(aa_gen[which(aa_gen[,4] == connector_labels[i]),5]), collapse = '')
  }
  binary_nss <- c()
  for(i in 1:length(connector_labels)){
    binary_nss[i] <- main_seq == connector_seqs[i]
  }
  
  connector_two_seq <- c()
  for(i in 1:length(binary_nss)){
    connector_two_seq[i] <- aa_gen[which(aa_gen[,4] == connector_labels[i]),5]
  }
  pos.135 <- c()
  for(i in 1:length(connector_two_seq)){
    pos.135[i] <- unlist(connector_two_seq[[i]])[135 + 16]
  }
    
}

##check 8's neighbors
{neighbors <- labels.new.new[which(dist_mat.large[,which(labels.new.new == labels.twenty[8])] == 2)]  
neighbors_seq <- c()
for(i in 1:length(neighbors)){
  neighbors_seq[i] <- aa_gen[which(aa_gen[,4] == neighbors[i]),5]
}
pos.131 <- c()
pos.142 <- c()
pos.261 <- c()
for(i in 1:length(twelve_neighbors)){
  pos.131[i] <- unlist(neighbors_seq[i])[131 + 16]
  pos.142[i] <- unlist(neighbors_seq[i])[142 + 16]
  pos.261[i] <- unlist(neighbors_seq[i])[261 + 16]
  
}
#####STOPPED HERE
x <- dist_mat.large[,which(neighbors[which(pos.142 == 'k')] == labels.new.new)]
main_connector_label <- labels.new.new[which(neighbors[which(pos.142 == 'k')] == labels.new.new)]
main_seq <- paste(unlist(aa_gen[which(aa_gen[,4] == main_connector_label),5]), collapse = '')
connector_labels <- labels.new.new[which(x == 1)]
connector_seqs <- c()
for(i in 1:length(connector_labels)){
  connector_seqs[i] <- paste(unlist(aa_gen[which(aa_gen[,4] == connector_labels[i]),5]), collapse = '')
}
binary_nss <- c()
for(i in 1:length(connector_labels)){
  binary_nss[i] <- main_seq == connector_seqs[i]
}

connector_two_seq <- c()
for(i in 1:length(binary_nss)){
  connector_two_seq[i] <- aa_gen[which(aa_gen[,4] == connector_labels[i]),5]
}
pos.261 <- c()
for(i in 1:length(connector_two_seq)){
  pos.261[i] <- unlist(connector_two_seq[[i]])[261 + 16]
}
}

#make an adj_mat from the dist mat large
date_indexes <- c()
for(i in 1:length(labels.new.new)){
  date_indexes[i] <- which(dna_gen[,4] == labels.new.new[i])
}
dates <- c()
for(i in 1:length(date_indexes)){
  index <- date_indexes[i]
  dates[i] <- paste(dna_gen[i,1], '-', dna_gen[i,2], '-', dna_gen[i,3], sep = '')
}


all_seqs <- c()
for(i in 1:length(labels.new.new)){
  all_seqs[i] <- aa_gen[which(aa_gen[,4] == labels.new.new[i]),5]
}
pos.121.all <- c()
for(i in 1:length(all_seqs)){
  pos.121.all[i] <- unlist(all_seqs[i])[121 + 16]
}
pos.121.all <- unlist(pos.121.all)

##get ratio for all 121k
{dist_mat_121k <- dist_mat.large[which(pos.121.all == 'k'), which(pos.121.all == 'k')]
labels_121k <- labels.new.new[which(pos.121.all == 'k')]
all_seqs_121k <- all_seqs[which(pos.121.all == 'k')]
dates_121k <- dates[which(pos.121.all == 'k')]
ones_121k <- list()
for(i in 1:ncol(dist_mat_121k)){
  ones_121k[[i]] <- labels_121k[which(dist_mat_121k[,i] == 1)]
}
lengths_121k <- c()
for(i in 1:length(ones_121k)){
 lengths_121k[i] <- length(unlist(ones_121k[[i]]))
}
remove <- which(lengths_121k < 5)
labels_121k <- labels_121k[-remove]
dist_mat_121k <- dist_mat_121k[-remove, -remove]
ones_121k <- ones_121k[-remove]
all_seqs_121k <- all_seqs_121k[-remove]
dates_121k <- dates_121k[-remove]
binary.121k <- list()
for(i in 1:length(ones_121k)){
  main <- paste(unlist(all_seqs_121k[i]), collapse = '')
  cluster <- ones_121k[[i]]
  binary_vec <- c()
  for(j in 1:length(cluster)){
    to_seq <- aa_gen[(which(cluster[j] == aa_gen[,4])),5]
    binary_vec[j] <- paste(unlist(to_seq), collapse = '') == main
  }
  binary.121k[[i]] <- binary_vec
}
ratios.121k <- c()
for(i in 1:length(binary)){
  ratios.121k[i] <- length(which(binary.121k[[i]] == FALSE)) / length(binary.121k[[i]])
}

}

##get ratio for all n121
{
  dist_mat_121n <- dist_mat.large[which(pos.121.all == 'n'), which(pos.121.all == 'n')]
  labels_121n <- labels.new.new[which(pos.121.all == 'n')]
  all_seqs_121n <- all_seqs[which(pos.121.all == 'n')]
  dates_121n <- dates[which(pos.121.all == 'n')]
  ones_121n <- list()
  for(i in 1:ncol(dist_mat_121n)){
    ones_121n[[i]] <- labels_121n[which(dist_mat_121n[,i] == 1)]
  }
  lengths_121n <- c()
  for(i in 1:length(ones_121n)){
    lengths_121n[i] <- length(unlist(ones_121n[[i]]))
  }
  remove <- which(lengths_121n < 5)
  labels_121n <- labels_121n[-remove]
  dist_mat_121n <- dist_mat_121n[-remove, -remove]
  ones_121n <- ones_121n[-remove]
  all_seqs_121n <- all_seqs_121n[-remove]
  dates_121n <- dates_121n[-remove]
  binary.121n <- list()
  for(i in 1:length(ones_121n)){
    main <- paste(unlist(all_seqs_121n[i]), collapse = '')
    cluster <- ones_121n[[i]]
    binary_vec <- c()
    for(j in 1:length(cluster)){
      to_seq <- aa_gen[(which(cluster[j] == aa_gen[,4])),5]
      binary_vec[j] <- paste(unlist(to_seq), collapse = '') == main
    }
    binary.121n[[i]] <- binary_vec
  }
  ratios.121n <- c()
  for(i in 1:length(binary.121n)){
    ratios.121n[i] <- length(which(binary.121n[[i]] == FALSE)) / length(binary.121n[[i]])
  }
  }

par(mfrow = c(1,1))
plot(as.Date(dates_121n, "%Y-%m-%d"), ratios.121n, pch = 16)
points(as.Date(dates_121k, "%Y-%m-%d"), ratios.121k, pch = 16, col = 'red')
legend('topright', legend = c('121n', '121k'), fill = c('black', 'red'))


dist_mat_k_test <- dist_mat.large[, which(pos.121.all == 'k')]
dist_mat_n_test <- dist_mat.large[,  which(pos.121.all == 'n')]
ones_k_test <- list()
lengths_k <- c()
ones_n_test <-list()
lengths_n <- c()
for(i in 1:ncol(dist_mat_k_test)){
  ones_k_test[[i]] <- which(dist_mat_k_test[,i] == 1)
  lengths_k[i] <- length(ones_k_test[[i]])
}
remove_k <- which(lengths_k < 5)
ones_k_test <- ones_k_test[-remove_k]

for(i in 1:ncol(dist_mat_n_test)){
  ones_n_test[[i]] <- which(dist_mat_n_test[,i] == 1)
  lengths_n[i] <- length(ones_n_test[[i]])
  
}
remove_n <- which(lengths_n < 5)
ones_n_test <- ones_n_test[-remove_n]

intersect(unlist(ones_k_test), unlist(ones_n_test))





