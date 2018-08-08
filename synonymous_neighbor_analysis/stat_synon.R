#####################
#Doc set up
#####################
setwd('~/Desktop/CSHL/')
load("most_prev_seq_by_year.RData")
load('AA_GEN.RData')
load('UNIQUE_DNA.RData')
library(lubridate)

#####################
#Set up names of amino acids
#####################
{{f <- c('ttt', 'ttc', rep(NA, 4))
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

seq <- c()
for(i in 1:length(prev)){
  seq[i] <- (aa_gen[which(unlist(aa_gen[,4]) == prev[i]), 5])
}
prev <- rbind(prev, seq)

synon_rate <- rep(0, 6)
for(i in 1:length(seq)){
  seq_this <- unlist(seq[i])
  for(j in 1:length(seq_this)){
    synon_rate[i] <- (unname(synon[which(names_aa == seq_this[j])]) + synon_rate[i]) 
  }
}
synon_rate <- synon_rate / (566 * 9)
non_synon_rate <- 1 - synon_rate
}

#####################
#Compare observed rates to theoretical, neighbor level 
#####################
{years <- c(2012, 2013, 2014, 2015, 2016, 2017)
list.neighbors <- list()
for(i in 1:length(years)){
  load(paste('network_', years[i], '.RData', sep = ''))
  assign(paste('network_', years[i], sep = ''), test)
  indexes <- get(paste('network_', years[i], sep = ''))[which(names(get(paste('network_', years[i], sep = ''))) == prev[1,i])]
  list.neighbors[[i]] <- names(get(paste('network_', years[i], sep = '')))[unlist(indexes)]
  remove(test)
}


binary.synon <- list()
for(i in 1:ncol(prev)){
  main <- which(aa_gen[,4] == unlist(prev[1,i]))
  list <- c()
  for(j in 1:length(unlist(list.neighbors[i]))){
    neighbor <- which(aa_gen[,4] == list.neighbors[[i]][j])
    list[j] <- paste(unlist(aa_gen[main,5]), collapse = '') == paste(unlist(aa_gen[neighbor,5]), collapse = '')
  }
  binary.synon[[i]] <- list
}

obs_synon <- c()
for(i in 1:ncol(prev)){
  obs_synon[i] <- length(which(binary.synon[[i]] == TRUE)) / length(binary.synon[[i]])
}

t.test(obs_synon, synon_rate)

}

#####################
#Compare observed rates to theoretical, annual level 
#####################
# {years <- c(2012, 2013, 2014, 2015, 2016, 2017)
# list.neighbors.year <- list()
# for(i in 1:length(years)){
#   load(paste('network_', years[i], '.RData', sep = ''))
#   list.neighbors.year[[i]] <- names(test)
#   remove(test)
# }
# 
# 
# 
# binary.synon.year <- list()
# for(i in 1:ncol(prev)){
#   main <- which(aa_gen[,4] == unlist(prev[1,i]))
#   list <- c()
#   for(j in 1:length(unlist(list.neighbors.year[i]))){
#     neighbor <- which(aa_gen[,4] == list.neighbors.year[[i]][j])
#     list[j] <- paste(unlist(aa_gen[main,5]), collapse = '') == paste(unlist(aa_gen[neighbor,5]), collapse = '')
#   }
#   binary.synon.year[[i]] <- list
#   print(i)
# }
# 
# obs_synon.year <- c()
# for(i in 1:ncol(prev)){
#   obs_synon.year[i] <- length(which(binary.synon.year[[i]] == TRUE)) / length(binary.synon.year[[i]])
# }
# 
# t.test(obs_synon.year, synon_rate)
# 
# }


networks <- c("network_2012", "network_2013", "network_2014", 
              "network_2015", "network_2016", "network_2017")
list.mass <- list()
for(j in 1:length(networks)){
  list.mass.network <- list()
  this <- get(networks[j])
  for(i in 1:length(this)){
    x <- this[which(names(this) == list.neighbors[[j]][i])]
    neighbors_neighbors <- names(this)[unlist(x)]
    list.mass.network[[i]] <- neighbors_neighbors
  }
  most_prev_neighbors <- names(this)[unlist(this[which(names(this) == prev[1,j])])]
  list.mass.network <- unlist(list.mass.network)
  list.mass[[j]] <- unique(c(list.mass.network, most_prev_neighbors))
}

##clean up unique

binary.synon.year <- list()
for(i in 1:ncol(prev)){
  main <- which(aa_gen[,4] == unlist(prev[1,i]))
  list <- c()
  for(j in 1:length(unlist(list.mass[i]))){
    neighbor <- which(aa_gen[,4] == list.mass[[i]][j])
    list[j] <- paste(unlist(aa_gen[main,5]), collapse = '') == paste(unlist(aa_gen[neighbor,5]), collapse = '')
  }
  binary.synon.year[[i]] <- list
  print(i)
}


obs_synon.year <- c()
for(i in 1:ncol(prev)){
  obs_synon.year[i] <- length(which(binary.synon.year[[i]] == TRUE)) / length(binary.synon.year[[i]])
}

t.test(obs_synon.year, synon_rate)

######FIND ALL NON SYNON MUTATIONS
mutations.synon.year <- list()
for(i in ncol(prev)){
  main <- which(aa_gen[,4] == unlist(prev[1,i]))
  main.seq <- unname(unlist(aa_gen[main,5]))
  temp.mut.list <- list()
  for(j in 1:length(unlist(list.mass[i]))){
    neighbor <- which(aa_gen[,4] == list.mass[[i]][j])
    neighbor.seq <- unname(unlist(aa_gen[neighbor,5]))
    mutations <- which(main.seq !=  neighbor.seq)
    if(length(mutations) == 0){next}
    mutations.vec <- c()
    for(k in 1:length(mutations)){
      mutations.vec[k] <- paste(main.seq[mutations[k]], mutations[k], neighbor.seq[mutations[k]], sep = ' ')
    }
    temp.mut.list[[j]] <- mutations.vec
  }
  mutations.synon.year[[i]] <- unique(unlist(temp.mut.list))
}
names(mutations.synon.year) <- prev[1,]

dates <- paste(aa_gen[,1], aa_gen[,2], aa_gen[,3], sep="-")
dates <- as.Date(dates, "%Y-%m-%d")
dates <- sort(dates)
dates <- unique(dates)


rug_2012 <- matrix(NA, nrow = length(dates), ncol = (length(unique(mutations.synon.year[[1]])) + 1))
rug_2012[,1] <- as.character(dates)
colnames(rug_2012) <- c('DATES', mutations.synon.year[[1]])
pos <- as.numeric(unlist(lapply(strsplit(mutations.synon.year[[1]], split = ' '), '[[', 2 )))
mut <- unlist(lapply(strsplit(mutations.synon.year[[1]], split = ' '), '[[', 3 ))
for(i in 1:length(mutations.synon.year[[1]])){
  for(j in 1:length(dates)){
    year <- year(dates[j])
    month <- month(dates[j])
    day <- day(dates[j])
    which.date.match <- which(aa_gen[,1] == year & aa_gen[,2] == month & aa_gen[,3] == day)
    seqs <- aa_gen[which.date.match,5]
    for(k in 1:length(seqs)){
      if(unlist(seqs[[k]])[pos[i]] == mut[i] | is.na(unlist(seqs[[k]])[pos[i]])){rug_2012[j,(i+1)] <- 1}else{rug_2012[j,(i+1)] <- 0}
      print(paste(i/length(mutations.synon.year[[1]]), j/length(dates), k/length(seqs), sep = '_'))
    }
  }
}
save(rug_2012, file = 'rug_2012.RData')

mut_names.2012 <- unlist(mutations.synon.year[[1]])
plot(x, y, pch = '|', xlim = c(0,1), ylim = c(0, length(mut_names.2012 / 2)),
     axes = FALSE, main = "Observed Point Mutations, 2012")
mtext(text = mut_names.2012, side = 2, line = 0.3, at = seq(1,length(mut_names.2012), length.out = length(mut_names.2012)), las = 2, cex=0.9)
mtext(text = years, side = 1, line = 0.3, at = seq(0,1, length.out = length(years)), las = 0, cex = 0.9)
box()

for(i in 4:dim(rug_2012)[2]){
  #rug_2012[,i] <- as.numeric(rug_2012[,i]) * i
  test <- (as.numeric(rug_2012[,i]) * i) - 1
  #plot.which <- which(rug_2012[,i] != 0)
  plot.which <- which(test != -1)
  x <- plot.which / length(dates)
  #y <- rug_2012[plot.which, i]
  y <- test[plot.which]
  points(x, y, pch = '|')
  print(i / dim(rug_2012)[2])
  
}

