setwd("~/Desktop/CSHL/")
library(seqinr)
library(fields)
library(colorRamps)

dms <- read.csv("a_perth_16_2009_H3N2.csv")
site <- dms[,1]
site <- as.character(site)
site_first <- as.integer(site[1:8])
site_last <- as.numeric(site[9:length(site)])


amino_acids <- c("K", "R", "H", "E", "D", "N", "Q", 
                 "T", "S", "C", "G", "A", "V", "L", 
                 "I", "M", "P", "Y", "F", "W")


fasta <- list()
fasta <- read.fasta("gisaid_epiflu_aa_2009.fasta")
#fasta <- read.fasta("gisaid_epiflu_aa_2007_2011.fasta")
fasta.ls <- list()
#list of lists of sequences
fasta.final <- list()
#vector of corresponding dates
fasta.names <- vector()
for(i in 1:length(fasta)){
  fasta.ls[[i]] <- unlist(fasta[i])
  fasta.final[[i]] <- unname(fasta.ls[[i]])
  #fasta.final[[i]] <- paste(fasta.final[[i]], collapse = '')
  # name <- names(fasta.ls[[i]])[1]
  # test <- strsplit(name, "")
  # test <- unlist(test)
  # test <- test[1:(length(test) -1)]
  # name <- paste(test, collapse = '')
  # fasta.names[i] <- name
}

#PROBLEM: SOME X's ARE LONGER THAN OTHERS
lengths <- rep(NA, length(fasta))
for(i in 1:length(fasta)){
  lengths[i] <- length(fasta[[i]])
}
mat.aa <- matrix(0, nrow = max(lengths), ncol = length(amino_acids))
colnames(mat.aa) <- amino_acids
for(i in 1:length(fasta)){
  x <- fasta.final[[i]]
  for(j in 1:max(lengths)){
    tryCatch({
      #A
      ifelse(x[j] == "k", mat.aa[j,1] <- (mat.aa[j,1] + 1), mat.aa[j,1])
      #G
      ifelse(x[j] == "r", mat.aa[j,2] <- (mat.aa[j,2] + 1), mat.aa[j,2])
      #C
      ifelse(x[j] == "h", mat.aa[j,3] <- (mat.aa[j,3] + 1), mat.aa[j,3])
      #T
      ifelse(x[j] == "e", mat.aa[j,4] <- (mat.aa[j,4] + 1), mat.aa[j,4])
      #T
      ifelse(x[j] == "d", mat.aa[j,5] <- (mat.aa[j,5] + 1), mat.aa[j,5])
      #T
      ifelse(x[j] == "n", mat.aa[j,6] <- (mat.aa[j,6] + 1), mat.aa[j,6])
      #T
      ifelse(x[j] == "q", mat.aa[j,7] <- (mat.aa[j,7] + 1), mat.aa[j,7])
      #T
      ifelse(x[j] == "t", mat.aa[j,8] <- (mat.aa[j,8] + 1), mat.aa[j,8])
      #T
      ifelse(x[j] == "s", mat.aa[j,9] <- (mat.aa[j,9] + 1), mat.aa[j,9])
      #T
      ifelse(x[j] == "c", mat.aa[j,10] <- (mat.aa[j,10] + 1), mat.aa[j,10])
      #T
      ifelse(x[j] == "g", mat.aa[j,11] <- (mat.aa[j,11] + 1), mat.aa[j,11])
      #T
      ifelse(x[j] == "a", mat.aa[j,12] <- (mat.aa[j,12] + 1), mat.aa[j,12])
      #T
      ifelse(x[j] == "v", mat.aa[j,13] <- (mat.aa[j,13] + 1), mat.aa[j,13])
      #T
      ifelse(x[j] == "l", mat.aa[j,14] <- (mat.aa[j,14] + 1), mat.aa[j,14])
      #T
      ifelse(x[j] == "i", mat.aa[j,15] <- (mat.aa[j,15] + 1), mat.aa[j,15])
      #T
      ifelse(x[j] == "m", mat.aa[j,16] <- (mat.aa[j,16] + 1), mat.aa[j,16])
      #T
      ifelse(x[j] == "p", mat.aa[j,17] <- (mat.aa[j,17] + 1), mat.aa[j,17])
      #T
      ifelse(x[j] == "y", mat.aa[j,18] <- (mat.aa[j,18] + 1), mat.aa[j,18])
      #T
      ifelse(x[j] == "f", mat.aa[j,19] <- (mat.aa[j,19] + 1), mat.aa[j,19])
      #T
      ifelse(x[j] == "w", mat.aa[j,20] <- (mat.aa[j,20] + 1), mat.aa[j,20])
      print(i)
    }, error=function(e){})
  }
  print(i)
}
mat.aa.prop <- matrix(0, nrow = max(lengths), ncol = length(amino_acids))
colnames(mat.aa.prop) <- amino_acids
for(i in 1:dim(mat.aa)[1]){
  total <- sum(mat.aa[i,])
  mat.aa.prop[i,] <- mat.aa[i,] / total
}

{par(mfrow = c(1,1), xpd = TRUE)
  image(mat.aa.prop,  col = c('lightgrey',matlab.like(565)),
        xlab = "Amino Acid Position",
        ylab = "Amino Acid",
        axes = F)
  title("Frequency of Each Amino Acid,\nH3 2009")
  mtext(text = amino_acids, side=2, line=0.3, at = seq(0,1, length.out = 20),  las=1, cex=0.8)
  mtext(text=seq(-8, 328, by = 20), side=1, line=0.3, at = seq(0,0.6, length.out = 17), las=2, cex=0.8, col = "DARKGREEN")
  mtext(text=seq(1, 222, by = 20), side=1, line=0.3, at = seq(0.65,1, length.out = 12), las=2, cex=0.8, col = "DARKBLUE")
  
  image.plot(mat.aa.prop, main = "Frequency of Each Amino Acid 
             /nH3 2008",
             xlab = "Amino Acid Position",
             ylab = "Amino Acid", legend.only = T)}
#abline(v = 0.7348152, col = "white", lwd = 2)

#######
#A/Perth/16/2009 H3N2 Data
#######
dms <- dms[,2:ncol(dms)]
dms <- as.matrix(dms)

dms_rearranged <- matrix(NA, nrow = dim(dms)[1], ncol = dim(dms)[2])
dms_rearranged[,1] <- dms[,9]
dms_rearranged[,2] <- dms[,15]
dms_rearranged[,3] <- dms[,7]
dms_rearranged[,4] <- dms[,4]
dms_rearranged[,5] <- dms[,3]
dms_rearranged[,6] <- dms[,12]
dms_rearranged[,7] <- dms[,14]
dms_rearranged[,8] <- dms[,17]
dms_rearranged[,9] <- dms[,16]
dms_rearranged[,10] <- dms[,2]
dms_rearranged[,11] <- dms[,6]
dms_rearranged[,12] <- dms[,1]
dms_rearranged[,13] <- dms[,18]
dms_rearranged[,14] <- dms[,10]
dms_rearranged[,15] <- dms[,8]
dms_rearranged[,16] <- dms[,11]
dms_rearranged[,17] <- dms[,13]
dms_rearranged[,18] <- dms[,20]
dms_rearranged[,19] <- dms[,5]
dms_rearranged[,20] <- dms[,19]


image(dms_rearranged,  col = matlab.like(566),
      xlab = "Amino Acid Position",
      ylab = "Amino Acid",
      axes = F)
title("DMS Frequency of Each Amino Acid,\nA/Perth/16/2009 (H3N2)")
mtext(text = amino_acids, side=2, line=0.3, at = seq(0,1, length.out = 20),  las=1, cex=0.8)
mtext(text=seq(-8, 328, by = 20), side=1, line=0.3, at = seq(0,0.6, length.out = 17), las=2, cex=0.8, col = "DARKGREEN")
mtext(text=seq(1, 222, by = 20), side=1, line=0.3, at = seq(0.65,1, length.out = 12), las=2, cex=0.8, col = "DARKBLUE")
image.plot(dms_rearranged, main = "Frequency of Each Amino Acid 
           /nH3 2008", 
           xlab = "Amino Acid Position",
           ylab = "Amino Acid", legend.only = T)

