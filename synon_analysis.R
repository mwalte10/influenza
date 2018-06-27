######################################
#BUILD SOME KIND OF GRAPHIC
#####################################
dna_list <- list(dna_2012_seq, dna_2013_seq, dna_2014_seq,
                 dna_2015_seq, dna_2016_seq, dna_2017_seq)
aa_list <- list(aa_2012_seq, aa_2013_seq, aa_2014_seq,
                aa_2015_seq, aa_2016_seq, aa_2017_seq)

dna_dif <- c()
aa_dif <- c()
for(i in 1:length(dna_list)){
  dna_dif[i] <- length(which(dna_list[[i]] != dna_list[[i + 1]]))
  aa_dif[i] <- length(which(aa_list[[i]] != aa_list[[i + 1]]))
}

total_synon <- (1 * 2) + (2 * 9) + (3 * 1) + (4 * 4) +  (6 * 3) 

dna_gen[,1] == 2012
length(which(dna_gen[,1] == 2014 & dna_gen[,7] == paste(aa_2014_seq, collapse = '')))

length(which(dna_gen[,7] == paste(aa_2012_seq, collapse = '')))


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