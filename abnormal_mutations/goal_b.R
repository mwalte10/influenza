setwd('~/Desktop/CSHL/')
dms_data <- read.csv('dms_data_prefs.csv')
mutations <- get(load('mutation_list.RData'))
load('AA_GEN.RData')
load('DNA_GEN.RData')
library(lubridate)
years <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)



abnormal_mutations <- list()
for(i in 1:length(mutations)){
  mut <- mutations[[i]]
  vec <- c()
  for(j in 1:length(mut)){
    site <- unlist(strsplit(mut[j], split = " "))[2]
    mutation <- unlist(strsplit(mut[j], split = " "))[3]
    if(is.na(match(toupper(mutation), colnames(dms_data)))){vec[j] <- NA && next}
    col <- which(colnames(dms_data) == toupper(mutation))
    ifelse(dms_data[site, col] < 0.01, vec[j] <- "ABNORMAL", vec[j]<- "NORMAL")
  }
  abnormal_mutations[[i]] <- vec
}
abnormal_mutations[[6]] <- c(abnormal_mutations[[6]], rep(NA, 2))

abnormal <- unlist(mutations)[which(unlist(abnormal_mutations) == "ABNORMAL")]
abnormal <- unname(abnormal)
positions <- c()
mutation <- c()
normal <- c()
for(i in 1:length(abnormal)){
	positions[i] <- strsplit(abnormal, split = ' ')[[i]][2]
	mutation[i] <- strsplit(abnormal, split = ' ')[[i]][3]
	normal[i] <- strsplit(abnormal, split = ' ')[[i]][1]
}
positions <- as.integer(positions)

normal <- normal[c(16, 37 ,13, 24)]
positions <- positions[c(16, 37 ,  13, 24)]
mutation <- mutation[c(16, 37 ,  13, 24)]


###############Create new Present/not present
abnormal_list.presence <- list()
dates <- paste(aa_gen[,1], aa_gen[,2], aa_gen[,3], sep="-")
dates <- as.Date(dates, "%Y-%m-%d")
dates <- sort(dates)

seqs.by.date <- list()
for(i in 1:length(unique(dates))){
	match.date <- which(dates == unique(dates)[i])
	seqs <- list()
	for(j in 1:length(match.date)){
		seqs[[j]] <- unname(unlist(aa_gen[match.date[j],5]))
	}
	seqs.by.date[[i]] <- seqs
	print(i / length(unique(dates)))
}


for(m in 1:length(positions)){
	vec <- c()
for(i in 1:length(seqs.by.date)){
	options <- seqs.by.date[[i]]
	test <- c()
	for(j in 1:length(options)){
	test[j] <- options[[j]][positions[m]]
	}
	test <- test[which(!is.na(test) == TRUE)]
	if(any(test == mutation[m])){vec[i] <- 1}else{vec[i] <- 0}
}
abnormal_list.presence[[m]] <- vec
print(m / length(abnormal))
}

##convert 1's to dates
date.list <- list()
for(i in 1:length(positions)){
	date.list[[i]] <- unique(dates)[which(abnormal_list.presence[[i]] == 1)]
}
save(date.list, file = "rug_dates_abnormal_muts.RData")


unique.mutation <- unique(normal)
codon.list <- list(c('gat', 'gac'),
c('aaa', 'aag'),
c('act', 'acc', 'aca'),
c('aca', 'acg'),
c('aaa', 'aag'),
c('act', 'acc', 'aca'),
c('aga', 'agg', 'cgt', 'cgc', 'cga', 'cgg'),
c('gat', 'gac'),
c('aca', 'acg'),
c('agt', 'agc'),
c('ctt', 'ctc', 'cta', 'ctg'),
c('aat', 'aac'),
c('ata'),
c('aat', 'aac'),
c('aat', 'aac'),
c('agt', 'agc'),
c('ctt', 'ctc'),
c('act', 'acc', 'aca', 'acg'),
c('aat', 'aac'),
c('gga', 'ggg'),
c('gat', 'gac'),
c('gat', 'gac'),
c('ggt', 'ggc'),
c('gta', 'gtg'),
c('act', 'acc'),
c('ttt', 'ttc'),
c('gat', 'gac'),
c('agt', 'agc'),
c('aaa', 'aag'),
c('aat', 'aac'),
c('gat', 'gac'),
c('aat', 'aac'),
c('aga', 'agg'),
c('tca', 'tcg'),
c('tat', 'tac'),
c('aat', 'aac'),
c('ggt', 'ggc'),
c('gaa', 'gag'))
codon.list <- codon.list[c(16, 37, 13, 24)]


dates.dna <- paste(dna_gen[,1], dna_gen[,2], dna_gen[,3], sep="-")
dates.dna <- as.Date(dates.dna, "%Y-%m-%d")
dates.dna <- sort(dates.dna)

available.date <- list()
for(i in 1:length(unique(dates.dna))){
	match.date <- which(dates.dna == unique(dates.dna)[i])
	seqs <- list()
	for(j in 1:length(match.date)){
		seqs[[j]] <- unname(unlist(dna_gen[match.date[j],6]))
	}
	available.date[[i]] <- seqs
	print(i / length(unique(dates.dna)))
}
#this isn't working
available <- list()
for(m in 1:length(positions)){
	vec <- c()
for(i in 1:length(available.date)){
	options <- available.date[[i]]
	test <- c()
	for(j in 1:length(options)){
	test[j] <- options[[j]][positions[m]]
	}
	test <- test[which(!is.na(test) == TRUE)]
	if(any(test == codon.list[[m]])){vec[i] <- 1}else{vec[i] <- 0}
}
available[[m]] <- vec
print(m / length(abnormal))
}
available.dates <- list()
for(i in 1:length(available)){
	available.dates[[i]] <- unique(dates.dna)[which(available[[i]] == 1)]
}
not_available <- list()
for(i in 1:length(available)){
	not_available[[i]] <- unique(dates.dna)[which(available[[i]] == 0)]
}

rug_plots_abnormal <- list(date.list, available.dates, not_available)
save(rug_plots_abnormal, file = "rug_plots_abnormal.RData")

###testing plotting
load('rug_plots_abnormal.RData')
date.list <- rug_plots_abnormal[1]
available.dates <- rug_plots_abnormal[2]
not_available <- rug_plots_abnormal[3]

spacing <- seq(0.4, 1.6, length.out = length(positions))
positions.lab <- c(124, 275, 140, 213)

ylabs <- paste(normal, positions.lab, mutation, sep = '')

{
par(mfrow = c(1,1), xpd = TRUE)
  plot(available.dates[[1]][[1]], y = rep(spacing[1], length(available.dates[[1]][[1]])), 
      pch = '|', col = "darkseagreen1", ylim = c(0.4, 1.6), ylab = '', xlab = '', yaxt='n')
  # plot(c(available.dates[[1]][[1]], not_available[[1]][[1]], date.list[[1]][[1]]), y = rep(0.5, length(c(available.dates[[1]][[1]], not_available[[1]][[1]], date.list[[1]][[1]]))), 
  #      pch = '|', 
  #      col = c(rep("darkseagreen1", length(available.dates[[1]][[1]])), rep('red', length(not_available[[1]][[1]])), rep('black', length(date.list[[1]][[1]]))), 
  #      ylim = c(0, 1), ylab = '', xlab = '', yaxt='n')
  
mtext(ylabs, side = 2, at = spacing, las = 2)
title(main = "Abnormal Substitutions Availability")
points(not_available[[1]][[1]], y = rep(spacing[1], length(not_available[[1]][[1]])), pch = '|', col = 'red')
points(date.list[[1]][[1]], y = rep(spacing[1], length(date.list[[1]][[1]])), pch = '|')
for(i in 2:length(spacing)){
	points(available.dates[[1]][[i]], y = rep(spacing[i], length(available.dates[[1]][[i]])), pch = '|', col = 'darkseagreen1')
	points(not_available[[1]][[i]], y = rep(spacing[i], length(not_available[[1]][[i]])), pch = '|', col = 'red')
	points(date.list[[1]][[i]], y = rep(spacing[i], length(date.list[[1]][[i]])), pch = '|')
}
legend(x = 13958.53, y = 1, legend = c('Not Accessible', "Observed", 'Accessbile'), 
       fill = c('red', 'black', 'darkseagreen1'), xpd = TRUE, cex = 0.6)
}

##prop of sequences with r142g

  vec <- c()
  for(i in 1:length(seqs.by.date)){
    options <- seqs.by.date[[i]]
    test <- c()
    for(j in 1:length(options)){
      test[j] <- options[[j]][142 + 16]
    }
    number.g <- length(which(test == 'g'))
    length.total <- length(test)
    vec[i] <- number.g / length.total
  }







for(i in 1:length(years)){
  load(paste('rug_', years[i], '.RData', sep = ''))
  assign(paste('rug_', years[i], sep = ''), rug)
  remove(rug)
}

mass_rug <- cbind(rug_2012, rug_2013, rug_2014,
                  rug_2015, rug_2016, rug_2017)
mass_rug_simp <- mass_rug[,c(1, match(abnormal, colnames(mass_rug)))]

dates <- paste(aa_gen[,1], aa_gen[,2], aa_gen[,3], sep="-")
dates <- as.Date(dates, "%Y-%m-%d")
dates <- sort(dates)
dates <- unique(dates)
mut_names <- colnames(mass_rug_simp)[2:ncol(mass_rug_simp)]

mass_rug_simp <- mass_rug_simp[,-5]
mut_names <- colnames(mass_rug_simp)[2:ncol(mass_rug_simp)]

i <- 2
test <- (as.numeric(mass_rug_simp[,i]) * i) - 1
plot.which <- which(test != -1)
x <- plot.which / length(dates)
y <- test[plot.which] 
plot(x, y, pch = '|', xlim = c(0,1), ylim = c(0, length(mut_names)),
     axes = FALSE, main = "Observed Point Mutations, \nAppears Rare from DMS")
mtext(text = mut_names, side = 2, line = 0.3, at = seq(1,length(mut_names), length.out = length(mut_names)), las = 2, cex=0.9)
mtext(text = years, side = 1, line = 0.3, at = seq(0,1, length.out = length(years)), las = 0, cex = 0.8)
box()
for(i in 3:dim(mass_rug_simp)[2]){
  test <- (as.numeric(mass_rug_simp[,i]) * i) - 1
  plot.which <- which(test != -1)
  x <- plot.which / length(dates)
  y <- test[plot.which] 
  points(x, y, pch = '|')
  print(i / dim(mass_rug_simp)[2])
}

which.point_291 <- c()
for(i in 1:dim(aa_gen)[1]){
  if(length(unlist(aa_gen[i, 5])) < 550){next}
  if(unname(unlist(aa_gen[i,5]))[291] == "d"){which.point_291[i] <- "MUT"}else{which.point_291[i] <- NA}
  print(i / nrow(aa_gen))
}
id_mut_291 <- unlist(aa_gen[which(which.point_291 == "MUT"),4])
id_mut_291.years <- c()
for(i in 1:length(id_mut_291)){
  id_mut_291.years[i] <- aa_gen[unlist(which(aa_gen[,4] == id_mut_291[i])), 1]
}
id_mut_291.years <- unlist(id_mut_291.years)
  

#########################
#Shows which ids/ years have the common mutations
#########################
which.point_156 <- c()
for(i in 1:dim(aa_gen)[1]){
  if(length(unlist(aa_gen[i, 5])) < 550){next}
  if(unname(unlist(aa_gen[i,5]))[156] == "r"){which.point_156[i] <- "MUT"}else{which.point_156[i] <- NA}
  print(i / nrow(aa_gen))
}  
id_mut_156 <- unlist(aa_gen[which(which.point_156 == "MUT"),4]) 
id_mut_156.years <- c()
for(i in 1:length(id_mut_156)){
  id_mut_156.years[i] <- aa_gen[unlist(which(aa_gen[,4] == id_mut_156[i])), 1]
}
id_mut_156.years <- unlist(id_mut_156.years)

id_mut <- unique(c(id_mut_291, id_mut_156))

  
  
  
