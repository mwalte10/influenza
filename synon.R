setwd('~/Desktop/CSHL/')

load('id.RData')
load('DNA.RData')
load('AA.RData')
load('date.RData')


dna_gen <- cbind(format(date, "%Y"),
                 format(date,"%m"),
                 format(date,"%d"),
                 id, dna_temp)
aa_gen <- cbind(format(date, "%Y"),
                format(date,"%m"),
                format(date,"%d"),
                id, aa_temp)


######################################
#SYNON MUTATIONS
#####################################
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
                     stop)}

dna_gen <- cbind(dna_gen, rep(NA, dim(dna_gen)[1]))
dna_test <- list()
for(i in 1:dim(dna_gen)[1]){
  x <- unname(unlist(dna_gen[i,5]))
  length <- ceiling((length(x) /3))
  codons <- c()
  for(j in 1:length){
    k <- j - 1
    codons[j] <- paste(x[(1:3) + 3 * k], collapse = '')
  }
  dna_test[[i]] <- codons
  print(i / dim(dna_gen)[1])
}
dna_gen[,6] <- dna_test

dna_gen <- cbind(dna_gen, rep(NA, dim(dna_gen)[1]))
test_aa <- list()
for(i in 1:dim(dna_gen)[1]){
  x <- unlist(dna_gen[i,6])
  aa <- c()
  for(j in 1:(length(x))){
    tryCatch({
      codon <- x[j]
      aa[j] <- colnames(amino_acids)[ceiling(which(amino_acids == codon) / 6)]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  codon_pos <- which(is.na(aa) == TRUE)
  given_aa <- unlist(aa_gen[i,5])
  aa[codon_pos] <- given_aa[codon_pos]
  test_aa[[i]] <- aa
  print(i / dim(dna_gen)[1])
}

for(i in 1:length(test_aa)){
  stop <- which(test_aa[[i]] == 'stop')[1]
  new <- test_aa[[i]]
  if(!is.na(stop)){test_aa[[i]] <- new[1:(stop - 1)]}
}
dna_gen[,7] <- test_aa

for(i in 1:length(test_aa)){
  dna_gen[i,7] <- paste(unlist(test_aa[[i]]), collapse = '')
}

save(dna_gen, file = 'DNA_GEN.RData')
save(aa_gen, file = 'AA_GEN.RData')

#unique DNA sequence table
which.unique.dna <- which(duplicated(dna_temp) == TRUE)
unique.date <- date[-which.unique.dna]
unique.year <- format(date, "%Y")[-which.unique.dna]
unique.month <- format(date,"%m")[-which.unique.dna]
unique.day <- format(date,"%d")[-which.unique.dna]
unique.id <- id[-which.unique.dna]
unique.dna <- dna_temp[-which.unique.dna]

dna <- cbind(unique.year, unique.month, unique.day,  
             unique.id, unique.dna, rep(0, length(unique.id)))
dna[,6] <- as.numeric(as.character(dna[,6]))
dna <- as.matrix(dna)
colnames(dna) <- c('YEAR', 'MONTH', 'DAY',
                   'ID', 'SEQUENCE', 'COUNT')

#unique AA sequence table
which.unique.aa <- which(duplicated(aa_temp) == TRUE)
unique.date <- date[-which.unique.aa]
unique.year <- format(date, "%Y")[-which.unique.aa]
unique.month <- format(date,"%m")[-which.unique.aa]
unique.day <- format(date,"%d")[-which.unique.aa]
unique.id <- id[-which.unique.aa]
unique.aa <- aa_temp[-which.unique.aa]

aa <- cbind(unique.year, unique.month, unique.day,  
            unique.id, unique.aa, rep(0, length(unique.aa)))
colnames(aa) <- c('YEAR', 'MONTH', 'DAY',
                  'ID', 'SEQUENCE', 'COUNT')


######################################
#SEPERATE BY YEAR
######################################
dna_test <- c()
for(i in 1:length(dna_temp)){
  dna_test[i] <- paste(dna_temp[[i]], collapse = '')
}

dna_unique <- c()
unique <- unique(dna_temp)
for(i in 1:length(unique)){
  dna_unique[i] <- paste(unique[[i]], collapse = '')
}

for(i in 1:dim(dna)[1]){
  dna[i,6] <- length(which(dna_test == dna_unique[i]))
}

aa_test <- c()
for(i in 1:length(aa_temp)){
  aa_test[i] <- paste(aa_temp[[i]], collapse = '')
}

aa_unique <- c()
unique_aa <- unique(aa_temp)
for(i in 1:length(unique_aa)){
  aa_unique[i] <- paste(unique_aa[[i]], collapse = '')
}

for(i in 1:dim(aa)[1]){
  aa[i,6] <- length(which(aa_test == aa_unique[i]))
}


save(aa, file = 'UNIQUE_AA.RData')
save(dna, file = 'UNIQUE_DNA.RData')


######################################
#FIND MOST COMMON SEQ FOR EACH YEAR
#####################################
dna_2012 <- dna[(which(dna[,1] == 2012)),]
dna_2013 <- dna[(which(dna[,1] == 2013)),]
dna_2014 <- dna[(which(dna[,1] == 2014)),]
dna_2015 <- dna[(which(dna[,1] == 2015)),]
dna_2016 <- dna[(which(dna[,1] == 2016)),]
dna_2017 <- dna[(which(dna[,1] == 2017)),]


aa_2012 <- aa[(which(aa[,1] == 2012)),]
aa_2013 <- aa[(which(aa[,1] == 2013)),]
aa_2014 <- aa[(which(aa[,1] == 2014)),]
aa_2015 <- aa[(which(aa[,1] == 2015)),]
aa_2016 <- aa[(which(aa[,1] == 2016)),]
aa_2017 <- aa[(which(aa[,1] == 2017)),]

dna_2012_seq <- as.character(unlist(dna_2012[which.max(dna_2012[,6]), 5]))
dna_2013_seq <- as.character(unlist(dna_2013[which.max(dna_2013[,6]), 5]))
dna_2014_seq <- as.character(unlist(dna_2014[which.max(dna_2014[,6]), 5]))
dna_2015_seq <- as.character(unlist(dna_2015[which.max(dna_2015[,6]), 5]))
dna_2016_seq <- as.character(unlist(dna_2016[which.max(dna_2016[,6]), 5]))
dna_2017_seq <- as.character(unlist(dna_2017[which.max(dna_2017[,6]), 5]))
years <- c(2012, 2013, 2014, 2015, 2016, 2017)
for(i in 1:length(years)){
  save(paste('dna_', years[i], '_seq', sep = ''), file = paste('DNA', years[i], '_seq.RData', sep = ''))
}
save(dna_2012_seq, file = 'DNA_2012_seq.RData')
save(dna_2013_seq, file = 'DNA_2013_seq.RData')
save(dna_2014_seq, file = 'DNA_2014_seq.RData')
save(dna_2015_seq, file = 'DNA_2015_seq.RData')
save(dna_2016_seq, file = 'DNA_2016_seq.RData')
save(dna_2017_seq, file = 'DNA_2017_seq.RData')





aa_2012_seq <- as.character(unlist(aa_2012[which.max(aa_2012[,6]), 5]))
aa_2013_seq <- as.character(unlist(aa_2013[which.max(aa_2013[,6]), 5]))
aa_2014_seq <- as.character(unlist(aa_2014[which.max(aa_2014[,6]), 5]))
aa_2015_seq <- as.character(unlist(aa_2015[which.max(aa_2015[,6]), 5]))
aa_2016_seq <- as.character(unlist(aa_2016[which.max(aa_2016[,6]), 5]))
aa_2017_seq <- as.character(unlist(aa_2017[which.max(aa_2017[,6]), 5]))
save(aa_2012_seq, file = 'AA_2012_seq.RData')
save(aa_2013_seq, file = 'AA_2013_seq.RData')
save(aa_2014_seq, file = 'AA_2014_seq.RData')
save(aa_2015_seq, file = 'AA_2015_seq.RData')
save(aa_2016_seq, file = 'AA_2016_seq.RData')
save(aa_2017_seq, file = 'AA_2017_seq.RData')





