setwd('~/Desktop/CSHL/')
library('seqinr')
library('stringr')


dna_2012_2015 <- read.fasta('dna_2012_2015.fasta')
dna_2016 <- read.fasta('dna_2016.fasta')
dna_2017 <- read.fasta('dna_2017.fasta')

dna <- c(dna_2012_2015, dna_2016, dna_2017)
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
date <- date[-16417]
dna_temp <- dna_temp[-remove_dna]
dna_temp <- dna_temp[-16417]
dna_id <- dna_id[-remove_dna]
dna_id <- dna_id[-16417]


#total number of unique nucleotide sequences 21937
length(unique(dna_temp))


aa_2012_2015 <- read.fasta('aa_2012_2015.fasta')
aa_2016 <- read.fasta('aa_2016.fasta')
aa_2017 <- read.fasta('aa_2017.fasta')

aa <- c(aa_2012_2015, aa_2016, aa_2017)
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

#total number of unique amino acid sequences is 10105
length(unique(aa_temp))

unique_aa_id <- unique(aa_id)

id <- unique_aa_id


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






