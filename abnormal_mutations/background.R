setwd('~/Desktop/CSHL/')
load('DNA_GEN.RData')
load('AA_GEN.RData')
library(lubridate)


dates <- paste(dna_gen[,1], dna_gen[,2], dna_gen[,3], sep="-")
dates <- as.Date(dates, "%Y-%m-%d")
dates <- sort(dates)

codons.291.final <- c()
codons.156.final <- c()
for(i in 1:length(dates)){
  year <- year(dates[i])
  month <- month(dates[i])
  day <- day(dates[i])
  date.which <- which(dna_gen[,1] == year & dna_gen[,2] == month & dna_gen[,3] == day)
  codon.date <- dna_gen[date.which,6]
  codon_291.vec <- c()
  codon_156.vec <- c()
  for(j in 1:length(date.which)){
    codon.date.select <- codon.date[[j]]
    codon_291.vec[j] <- codon.date.select[291]
    codon_156.vec[j] <- codon.date.select[156]
  }
  if(any(codon_291.vec == "ggt" | codon_291.vec == "ggc"| codon_291.vec == "gat"| codon_291.vec == "gac")){codons.291.final[i] <- 1}else{codons.291.final[i] <- 0}
  if(any(codon_156.vec == "ata" | codon_156.vec == "cgt" | codon_156.vec == "cgc" | codon_156.vec == "cga" | codon_156.vec == "cgg" | codon_156.vec == "aga" | codon_156.vec == "agg")){codons.156.final[i] <- 1}else{codons.156.final[i] <- 0}
  print(i / length(dates))
}


rug_dates.291 <- c()
rug_dates.156 <- c()
for(i in 1:length(dates)){
  if(codons.291.final[i] == 1){rug_dates.291[i] <- 1}else{rug_dates.291[i] <- 0}
  if(codons.156.final[i] == 1){rug_dates.156[i] <- 1}else{rug_dates.156[i] <- 0}
  print(i / length(dates))
}

seen.156 <- c()
seen.291 <- c()
for(i in 1:nrow(dna_gen)){
  sequence <- unname(unlist(aa_gen[which(aa_gen[,4] == unlist(dna_gen[i,4])),5]))
  if(length(sequence) < 500){next}
  if(is.na(sequence[156])){next}
  if(is.na(sequence[291])){next}
  if(sequence[156] == 'r'){seen.156[i] <- 1}else{seen.156[i] <- 0}
  if(sequence[291] == 'd'){seen.291[i] <- 1}else{seen.291[i] <- 0}
  print(i / nrow(dna_gen))
}
seen.156[which(is.na(seen.156) == "TRUE")] <- 0
seen.291[which(is.na(seen.291) == "TRUE")] <- 0



##RUN RUGS FOR TWO MUTATIONS

#call dna_gen instead of aa_gen, put in positioning 

years <- c(2012, 2013, 2014, 2015, 2016, 2017)
for(i in 1:length(years)){
  load(paste('rug_', years[i], '.RData', sep = ''))
  assign(paste('rug_', years[i], sep = ''), rug)
  remove(rug)
}

rug.156 <- dates[which(seen.156 == 1)]
rug.156.avail <- dates[which(codons.156.final == 1)]
rug.156.navail <- dates[which(codons.156.final != 1)]

rug.291 <- dates[which(seen.291 == 1)]
rug.291.avail <- dates[which(codons.291.final == 1)]
rug.291.navail <- dates[which(codons.291.final != 1)]


##times seen isn't plotting
plot(x = rug.156.avail, y = rep(1, length(rug.156.avail)), col = "forestgreen", pch = "|", 
     main = "Availability of i156r Mutation")
points(x = rug.156, y = rep(1, length(rug.156)), pch = "|")
points(x = rug.156.navail, y = rep(1, length(rug.156.navail)), col = "red", pch = "|")
legend("topright", c("Available", "Not Available", "Observed"), fill = c('forestgreen', 'red', 'black'))


plot(x = rug.291.avail, y = rep(1, length(rug.291.avail)), col = "forestgreen", pch = "|", 
     main = "Availability of g291d Mutation")
points(x = rug.291, y = rep(1, length(rug.291)), pch = "|")
points(x = rug.291.navail, y = rep(1, length(rug.291.navail)), col = "red", pch = "|")
legend("topright", c("Available", "Not Available", "Observed"), fill = c('forestgreen', 'red', 'black'))












