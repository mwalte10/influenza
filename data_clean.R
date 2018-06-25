setwd('~/Desktop/CSHL/')
library('seqinr')
library('stringr')


dna_2012_2015 <- read.fasta('dna_2012_2015.fasta')
dna_2016 <- read.fasta('dna_2016.fasta')
dna_2017 <- read.fasta('dna_2017.fasta')

dna <- c(dna_2012_2015, dna_2016, dna_2017)
dna_temp <- list()
dna_id <- c()
for(i in 1:length(dna)){
  x <- unlist(dna[[i]])
  dna_temp[[i]] <- paste(x, collapse = '')
  beg <- str_locate(dna_temp[[i]], "atg")[1]
  dna_temp[[i]] <- paste(x[beg:length(x)], sep = '')
  a <- unlist(strsplit(unlist(getAnnot(x)), ' '))
  dna_id[i] <- a[length(a)]
}
unique_dna_id <- unique(dna_id)

#remove any sequences with the same ID
unique_dna <- dna_id[-which(duplicated(dna_id) == TRUE)]

#total number of nucleotide sequences is 39572
dna_temp <- dna_temp[-which(duplicated(dna_id) == TRUE)]
#total number of unique nucleotide sequences 22105
unique_dna <- unique(dna_temp)
length(unique_dna)


aa_2012_2015 <- read.fasta('aa_2012_2015.fasta')
aa_2016 <- read.fasta('aa_2016.fasta')
aa_2017 <- read.fasta('aa_2017.fasta')

aa <- c(aa_2012_2015, aa_2016, aa_2017)
aa_temp <- list()
aa_id <- c()
for(i in 1:length(aa)){
  tryCatch({
    x <- unlist(aa[[i]])
    aa_temp[[i]] <- paste(x, collapse = '')
    beg <- str_locate(aa_temp[[i]], "m")[1]
    aa_temp[[i]] <- paste(x[beg:length(x)], sep = '')
    aa_id[i] <- attr(x,"name")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#remove any sequences with the same ID
unique_aa_id <- aa_id[-which(duplicated(aa_id) == TRUE)]
#total number of aa sequences is 39572
aa_temp <- aa_temp[-which(duplicated(aa_id) == TRUE)]
length(aa)
#total number of unique amino acid sequences is 10170
unique_aa <- unique(aa)
length(unique_aa)

id <- unique_aa_id


write.fasta(dna_temp, unique_dna_id,'DNA.fasta')
write.fasta(aa, unique_aa_id,'AA.fasta')
save(id, file = 'id.RData')
save(dna_temp, file = 'DNA.RData')
save(aa_temp, file = 'AA.RData')
