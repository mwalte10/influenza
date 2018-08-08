setwd('~/Desktop/CSHL/')
library('seqinr')
library('stringr')


# dna_2012_2015 <- read.fasta('dna_2012_2015.fasta')
# dna_2016 <- read.fasta('dna_2016.fasta')
# dna_2017 <- read.fasta('dna_2017.fasta')
dna_2011 <- read.fasta('dna_2011.fasta')

# dna <- c(dna_2012_2015, dna_2016, dna_2017)
dna <- dna_2011
dna_temp <- list()
dna_id <- c()
date <- c()
for(i in 1:length(dna)){
  x <- unlist(dna[[i]])
  dna_temp[[i]] <- paste(x, collapse = '')
  beg <- str_locate(dna_temp[[i]], "atg")[1]
  dna_temp[[i]] <- paste(x[beg:length(x)], sep = '')
  a <- unlist(strsplit(unlist(getAnnot(x)), ' '))
  dna_id[i] <- a[length(a)]
  date[i] <- a[1]
}

#clean date
date_temp <- c()
for(i in 1:length(date)){
  date_temp[i] <- unlist(strsplit(unlist(date[i]), '>'))[2]
}
date <- date_temp
date <- as.Date(date)

#remove duplicates and those w/o dates
duplicates <- which(duplicated(dna_id) == TRUE)
no_dates <- which(is.na(date) == TRUE)
remove_dna <- c(duplicates, no_dates)

#total number of nucleotide sequences is 39304
date <- date[-remove_dna]
# date <- date[-16417]
dna_temp <- dna_temp[-remove_dna]
# dna_temp <- dna_temp[-16417]
dna_id <- dna_id[-remove_dna]
# dna_id <- dna_id[-16417]

dna_temp_triplet <- list()
for(j in 1:length(dna_temp)){
  aa <- c()
  for(i in 1:(566)){
    triplet <- seq(1,3) + (i - 1) * (3)
    aa[i] <- paste(dna_temp[[j]][triplet], collapse = '')
  }
  dna_temp_triplet[[j]] <- aa
  print(j / length(dna_temp))
}


f <- c('ttt', 'ttc', rep(NA, 4))
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


amino_acids.names <- c('f', 'l', 's', 'y',
                 'c', 'w', 'p', 'h',
                 'q', 'r', 'i', 'm',
                 't', 'n', 'k', 'v',
                 'a', 'd', 'e', 'g')

inferred_aa <- list()
for(i in 1:length(dna_temp_triplet)){
  triplet <- dna_temp_triplet[[i]]
  aa <- c()
  for(j in 1:length(triplet)){
    x <- which(amino_acids == triplet[j])
    if(length(x) == 0){next}
    aa[j] <- colnames(amino_acids)[ceiling(x / 6)]
  }
  inferred_aa[[i]] <- aa
}

dna_2010.mat <- cbind(year(date), month(date), day(date), dna_id, dna_temp, dna_temp_triplet, inferred_aa)

##This is set up
dna_2011.mat <- cbind(year(date), month(date), day(date), dna_id, dna_temp, dna_temp_triplet, inferred_aa)
dna_gen.new <- rbind(dna_2010.mat, dna_2011, dna_gen)
dna_gen <- dna_gen.new

dna_gen_2010 <- dna_gen[1:1280,]
dna_gen_other <- dna_gen[1839:nrow(dna_gen),]
dna_gen.test <- rbind(dna_gen_2010, dna_2011.mat, dna_gen_other)
dna_gen <- dna_gen.test

save(dna_gen, file = "DNA_GEN.RData")

dna_gen[39305:nrow(dna_gen), 6] <- dna_temp_triplet
dna_gen[39305:nrow(dna_gen), 7] <- inferred_aa


#total number of unique nucleotide sequences 21937
length(unique(dna_temp))
which.unique <- match(unique(dna_temp), dna_temp)

load('DNA_GEN.RData')
load('UNIQUE_DNA.RData')

add_2018 <- cbind(year(date), month(date), day(date), dna_id, dna_temp, rep(NA, length(dna_id)), rep(NA, length(dna_id)))
new_dna_gen <- rbind(dna_gen, add_2018)
dna_gen <- new_dna_gen

add_2018.unique <- cbind(year(date)[which.unique], month(date)[which.unique], day(date)[which.unique], 
                         dna_id[which.unique], dna_temp[which.unique], 
                         rep(NA, length(which.unique)), rep(NA, length(which.unique)))
dna <- rbind(dna, add_2018.unique)

save(dna_gen, file = 'DNA_GEN.RData')
save(dna, file = 'UNIQUE_DNA.RData')























# aa_2012_2015 <- read.fasta('aa_2012_2015.fasta')
# aa_2016 <- read.fasta('aa_2016.fasta')
# aa_2017 <- read.fasta('aa_2017.fasta')
aa_2010 <- read.fasta('aa_2010.fasta')

# aa <- c(aa_2012_2015, aa_2016, aa_2017)
aa <- aa_2010
aa_temp <- list()
aa_id <- c()
aa_date <- c()
for(i in 1:length(aa)){
  tryCatch({
    x <- unlist(aa[[i]])
    aa_temp[[i]] <- paste(x, collapse = '')
    beg <- str_locate(aa_temp[[i]], "m")[1]
    aa_temp[[i]] <- paste(x[beg:length(x)], sep = '')
    aa_id[i] <- attr(x,"name")
    a <- unlist(strsplit(unlist(getAnnot(x)), ' '))
    aa_date[i] <- a[length(a)] 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#remove duplicates and those w/o dates
aa_date <- as.Date(aa_date)
unknown_date <- which(is.na(aa_date) == TRUE)
duplicate_id <- which(duplicated(aa_id) == TRUE)

remove_aa <- c(unknown_date, duplicate_id)

#total number of aa sequences is 39304
aa_temp <- aa_temp[-remove_aa]
aa_date <- aa_date[-remove_aa]
aa_id <- aa_id[-remove_aa]
load('AA_GEN.RData')

add_2010 <- cbind(year(aa_date), month(aa_date), day(aa_date), aa_id, aa_temp)
add_2011 <- cbind(year(aa_date), month(aa_date), day(aa_date), aa_id, aa_temp)
aa_gen.new <- rbind(add_2010, add_2011, aa_gen)
aa_gen <- aa_gen.new
save(aa_gen, file = "AA_GEN.RData")

#total number of unique amino acid sequences is 10105
length(unique(aa_temp))


for(i in (10106:nrow(aa))){
  aa[i,5] <- paste(unlist(aa[i, 5]), collapse = '')
}

add_2018.aa <- cbind(year(aa_date), month(aa_date), day(aa_date), 
                     aa_id, aa_temp)
aa_gen <- rbind(aa_gen, add_2018.aa)
save(aa_gen, file = 'AA_GEN.RData')

which.unique <- match(unique(aa_temp), aa_temp)
add_2018.aa.unique <- cbind(year(aa_date)[which.unique], month(aa_date)[which.unique], day(aa_date)[which.unique], 
                     aa_id[which.unique], aa_temp[which.unique], rep(NA, length(which.unique)))
aa <- rbind(aa, add_2018.aa.unique)
save(aa, file = 'UNIQUE_AA.RData')


unique_aa_id <- unique(aa_id)

id <- unique_aa_id

{
####ORDER EVERYTHING BY DATE
date <- date[order(date)]
dna_temp <- dna_temp[order(date)]
aa_temp <- aa_temp[order(date)]
id <- id[order(date)]


write.fasta(dna_temp, unique_dna_id,'DNA.fasta')
write.fasta(aa, unique_aa_id,'AA.fasta')
save(id, file = 'id.RData')
save(dna_temp, file = 'DNA.RData')
save(aa_temp, file = 'AA.RData')
save(date, file = 'date.RData')

#CHECK TO ENSURE THAT AAs and DNA are the same
# for(i in 1:length(date)){
#   if(date[i] != aa_date[i]){print(i) & break}
# }

# for(i in 1:length(aa_id)){
#   if(aa_id[i] != dna_id[i]){print(i) & break}
# }

}




