######################################
#BUILD SOME KIND OF GRAPHIC
#####################################
# dna_list <- list(dna_2012_seq, dna_2013_seq, dna_2014_seq,
#                  dna_2015_seq, dna_2016_seq, dna_2017_seq)
# aa_list <- list(aa_2012_seq, aa_2013_seq, aa_2014_seq,
#                 aa_2015_seq, aa_2016_seq, aa_2017_seq)
# 
# dna_dif <- c()
# aa_dif <- c()
# for(i in 1:length(dna_list)){
#   dna_dif[i] <- length(which(dna_list[[i]] != dna_list[[i + 1]]))
#   aa_dif[i] <- length(which(aa_list[[i]] != aa_list[[i + 1]]))
# }

total_synon <- (1 * 2) + (2 * 9) + (3 * 1) + (4 * 4) +  (6 * 3) 

# dna_gen[,1] == 2012
# length(which(dna_gen[,1] == 2014 & dna_gen[,7] == paste(aa_2014_seq, collapse = '')))
# 
# length(which(dna_gen[,7] == paste(aa_2012_seq, collapse = '')))


# dna_2012 <- dna[(which(dna[,1] == 2012)),]
# dna_2013 <- dna[(which(dna[,1] == 2013)),]
# dna_2014 <- dna[(which(dna[,1] == 2014)),]
# dna_2015 <- dna[(which(dna[,1] == 2015)),]
# dna_2016 <- dna[(which(dna[,1] == 2016)),]
# dna_2017 <- dna[(which(dna[,1] == 2017)),]
# 
# 
# aa_2012 <- aa[(which(aa[,1] == 2012)),]
# aa_2013 <- aa[(which(aa[,1] == 2013)),]
# aa_2014 <- aa[(which(aa[,1] == 2014)),]
# aa_2015 <- aa[(which(aa[,1] == 2015)),]
# aa_2016 <- aa[(which(aa[,1] == 2016)),]
# aa_2017 <- aa[(which(aa[,1] == 2017)),]
# 
# 


###############
#LOAD IN DATA
###############
setwd('~/Desktop/CSHL/')
library(Matrix)
years <- c(2012, 2013, 2014, 2015, 2016, 2017)
for(i in 1:length(years)){
  load(paste('AA_', years[i], '_seq.RData', sep = ''))
  load(paste('DNA_', years[i], '_seq.RData', sep = ''))
}
load('AA_GEN.RData')
load('DNA_GEN.RData')
load('UNIQUE_AA.RData')
load('UNIQUE_DNA.RData')

###############
#FIND HOW MANY DNA SEQs REP ONE AA
###############
rep <- c()
for(i in 1:dim(aa)[1]){
  x <- unlist(aa[i,5])
  rep[i] <- length(which(dna[,7] == paste(x, collapse = '')))
}

#those which differ 
differ <- which(rep != aa[,6])
differ_mat <- cbind(rep[differ], unlist(aa[differ,6]))
colnames(differ_mat) <- c('UNIQUE_SEQ','COUNT')


###############
#SYNON MUT COUNTS
###############
# {f <- c('ttt', 'ttc', rep(NA, 4))
#   l <- c('tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg')
#   s <- c('tct', 'tcc', 'tca', 'tcg', 'agt', 'agc')
#   y <- c('tat', 'tac', rep(NA, 4))
#   c <- c('tgt', 'tgc', rep(NA, 4))
#   w <- c('tgg', rep(NA, 5))
#   p <- c('cct', 'ccc', 'cca', 'ccg', rep(NA, 2))
#   h <- c('cat', 'cac', rep(NA, 4))
#   q <- c('caa', 'cag', rep(NA, 4))
#   r <- c('cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg')
#   i <- c('att', 'atc', 'ata', rep(NA, 3))
#   m <- c('atg', rep(NA, 5))
#   t <- c('act', 'acc', 'aca', 'acg', rep(NA, 2))
#   n <- c('aat', 'aac', rep(NA, 4))
#   k <- c('aaa', 'aag', rep(NA, 4))
#   v <- c('gtt', 'gtc', 'gta', 'gtg', rep(NA, 2))
#   a <- c('gct', 'gcc', 'gca', 'gcg', rep(NA, 2))
#   d <- c('gat', 'gac', rep(NA, 4))
#   e <- c('gaa', 'gag', rep(NA, 4))
#   g <- c('ggt', 'ggc', 'gga', 'ggg', rep(NA, 2))
#   stop <- c('tga', 'taa', 'tag', rep(NA, 3))
#   amino_acids <- cbind(f, l, s, y,
#                        c, w, p, h,
#                        q, r, i, m,
#                        t, n, k, v,
#                        a, d, e, g,
#                        stop)}


amino_acids <- c('f', 'l', 's', 'y',
                     'c', 'w', 'p', 'h',
                     'q', 'r', 'i', 'm',
                     't', 'n', 'k', 'v',
                     'a', 'd', 'e', 'g')
amino_acids_count <- c(2, 6, 6, 2, 2, 1, 4, 2, 2, 6, 3, 1, 4, 2, 2, 4, 4, 2, 2, 4)

#getting inf values, meaning the number of synonymous mutations is just too big
synon <- rep(1, dim(aa)[1])
for(i in 1:dim(aa)[1]){
  aa_test <- unlist(aa[i,5])
  for(i in 1:length(aa_test)){
    count <- amino_acids_count[which(amino_acids == aa_test[j])]
    synon[i] <- synon[i] * count
  }
  print(i / dim(aa)[1])
}

aa_2012 <- aa[(which(aa[,1] == 2012)),]
aa_2013 <- aa[(which(aa[,1] == 2013)),]
aa_2014 <- aa[(which(aa[,1] == 2014)),]
aa_2015 <- aa[(which(aa[,1] == 2015)),]
aa_2016 <- aa[(which(aa[,1] == 2016)),]
aa_2017 <- aa[(which(aa[,1] == 2017)),]

dna_2012 <- dna[(which(dna[,1] == 2012)),]
dna_2013 <- dna[(which(dna[,1] == 2013)),]
dna_2014 <- dna[(which(dna[,1] == 2014)),]
dna_2015 <- dna[(which(dna[,1] == 2015)),]
dna_2016 <- dna[(which(dna[,1] == 2016)),]
dna_2017 <- dna[(which(dna[,1] == 2017)),]


###This is abandoned
{which.list <- list()
synon_rate <- c()
for(i in 1:length(years)){
  dna.year <- get(paste('dna_', years[i], sep = ''))
  aa.year <- get(paste('aa_', years[i], '_seq', sep = ''))
  which.list[[i]] <- which(dna.year[,7] == paste(aa.year, collapse = ''))
  aa.synon.obs <- list()
  for(j in 1:length(which.list[[i]])){
    z <- unlist(dna_2012[which.2012[j], 5])
    aa.synon.obs[[j]] <- ceiling(which(z != dna_2012_seq) / 3)
    print(j / length(which.2012))
  }
  aa.synon.obs.unique <- unique(unlist(aa.synon.obs))
  test <- aa_2012_seq[aa.synon.obs.unique]
  no.na <- test[!is.na(aa_2012_seq[aa.synon.obs.unique]) == TRUE]
  synon_2012 <- 1
  for(k in 1:length(no.na)){
    count <- amino_acids_count[which(amino_acids == no.na[k])]
    synon_year <- synon_year * count 
    synon_rate[i] <- synon_year 
  }
}

which.2012 <- which(dna_2012[,7] == paste(aa_2012_seq, collapse = ''))

aa.synon.obs <- list()
for(i in 1:length(which.2012)){
  z <- unlist(dna_2012[which.2012[i], 5])
  aa.synon.obs[[i]] <- ceiling(which(z != dna_2012_seq) / 3)
  print(i / length(which.2012))
}
#positions which see synonmous mutations
aa.synon.obs.unique <- unique(unlist(aa.synon.obs))
test <- aa_2012_seq[aa.synon.obs.unique]
no.na <- test[!is.na(aa_2012_seq[aa.synon.obs.unique]) == TRUE]
synon_2012 <- 1
for(i in 1:length(no.na)){
  count <- amino_acids_count[which(amino_acids == no.na[i])]
  synon_2012 <- synon_2012 * count 
}
synon_rate_2012 <- length(unlist(aa.synon.obs)) / synon_2012
}

adj_mat_2013 <- Matrix(NA, nrow = dim(dna_2013)[1], ncol = dim(dna_2013)[1], sparse = TRUE)
#save(adj_mat_2013, file = "adj_mat_2013.RData")

for(i in 1:dim(dna_2013)[1]){
  comp <- dna_2013[,5][i]
  comp <- unlist(comp)
  for(j in 1:dim(dna_2013)[1]){
    to <- dna_2013[,5][j]
    to <- unlist(to)
    x <- length(which(comp != to))
    if(x == 1){adj_mat_2013[i,j] <- 1} 
  }
  print(i / dim(dna_2013)[1])
}
colnames(adj_mat_2013) <- dna_2013[,4]
rownames(adj_mat_2013) <- dna_2013[,4]
id_2013 <- unlist(dna_2013[,4])

#find all sequences with neighbors
cols <- ceiling(which(adj_mat_2013 == 1) / 1671)
unique.cols <- unique(cols)
counts <- c()
#find how many neighbors each sequenec has
for(i in 1:length(unique.cols)){
  counts[i] <- length(which(cols == unique.cols[i]))
}
names(counts) <- id_2013[unique.cols]

rows <- (ceiling(which(adj_mat_2013 == 1) / 1671) - (which(adj_mat_2013 == 1) / 1671)) * 1671
rows[116] <- 250
new <- cbind(id_2013[cols], id_2013[rows])
test <- network(new, matrix.type = 'edgelist', directed = FALSE)
plot.network(test, main = "2013 DNA Sequence Network \nAll Unique Sequences", vertex.col = 'grey', 
             usecurve = 0, #vertex.cex = seq(0.2, 1.3, length.out = max(counts))
             jitter= FALSE)
#save(test, file ='test.RData')


###### NETWORK WITH SEQUENCES WITH MORE THAN ONE NEIGHBOR
which.not.one <- unname(which(counts != 1))
names_not_one <- names(counts[which.one])
keep_left <- list()
keep_right <- list()
for(i in 1:length(names_not_one)){
  keep_left[[i]] <- which(new[,1] == names_not_one[i])
  keep_right[[i]] <- which(new[,2] == names_not_one[i])
}
test_ah <- new[keep_left,]
test_ah <- rbind(test_ah, new[keep_right,])
test_sh <- network(test_ah, matrix.type = 'edgelist', directed = FALSE)
plot.network(test_sh, main = "2013 DNA Sequence Network \nMore Than One Neighbor", vertex.col = 'grey', 
             usecurve = 0, #vertex.cex = seq(0.2, 1.3, length.out = max(counts))
             jitter= FALSE)
#save(test_sh, file = 'test_more_than_one_neighbor.RData')


