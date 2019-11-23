# The following code is a parser for Max Quant results to be fed into Percolator
# and analysis of the results of percolator in the end.

#+------------------------------+----+
#| PREPARING THE EVIDENCE FILE: | RT |
#+------------------------------+----+

ms2 <- read.delim('C:/Users/nslavov/Desktop/New MSMS/msms.txt', header = TRUE, stringsAsFactors = FALSE)
ev <- read.delim('C:/Users/nslavov/Desktop/New MSMS/evidence.txt', header = TRUE, stringsAsFactors = FALSE)

# Adding features from MSMS file:
ev$Scan.index <- ms2$Scan.index[match(ev$Best.MS.MS, ms2$id)] #for matching duty cycles from MSMS Scans
ev$n.matches <- ms2$Number.of.Matches[match(ev$Best.MS.MS, ms2$id)]
ev$ScanNr <- ms2$Scan.number[match(ev$Best.MS.MS, ms2$id)]

write.table(ev, 'C:/Users/nslavov/Desktop/Percolator/MQ files/evidence+3 columns of MSMS.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

ev <- read.delim('C:/Users/nslavov/Desktop/Percolator/MQ files/evidence+3 columns of MSMS.txt', header = TRUE, stringsAsFactors = FALSE)

ev <- ev[, c("Raw.file", "Reverse", "ScanNr", "Retention.time", "Mass.Error..Da.", "PEP", "Length",
             "Mass", "m.z", "Score", "Delta.score", "Resolution", "Missed.cleavages", "n.matches",
             "Sequence","Leading.razor.protein", "Modified.sequence", "Charge", "Scan.index", "id")]

ev$Raw.file <- paste(ev$Raw.file, "_ID", ev$id, "_SI", ev$Scan.index, sep="")
ev$Raw.file <- paste(ev$Raw.file, sep="")
ev$id <- NULL
ev$Scan.index <- NULL
colnames(ev)[1] <- "id"

ev[!ev$Reverse=="+", "Reverse"] <- 1
ev[ev$Reverse=="+", "Reverse"] <- -1
colnames(ev)[2] <- "label"

colnames(ev)[4] <- "RT"
colnames(ev)[7] <- "PepLen"
colnames(ev)[5] <- "dM"
colnames(ev)[11] <- "dScore"
colnames(ev)[13] <- "enzInt"
colnames(ev)[15] <- "peptide"
colnames(ev)[16] <- "proteinId1"

ev$theoMass <- ev$Mass + ev$Charge*1.007276 # + 229.162932*ev$nTMT
ev$obsMass <- ev$m.z * ev$Charge
ev$dM <- abs(ev$theoMass-ev$obsMass)
ev$theoMass <- NULL
ev$obsMass <- NULL

ev$Charge1 <- 0
ev[ev$Charge==1, "Charge1"] <- 1
ev$Charge2 <- 0
ev[ev$Charge==2, "Charge2"] <- 1
ev$Charge3 <- 0
ev[ev$Charge==3, "Charge3"] <- 1
ev$Charge4 <- 0
ev[ev$Charge==4, "Charge4"] <- 1
ev$Charge5 <- 0
ev[ev$Charge==5, "Charge5"] <- 1
ev$Charge6 <- 0
ev[ev$Charge==6, "Charge6"] <- 1
ev$Charge <- NULL

# Separating Reverse vs. PSMs vs. Contaminent:
# Many of the reverse sequences are not marked as reversed:
#ev.rev <- ev[ev$label==-1, ]  # >>> 281,302
#ev.psm <- ev[ev$label==1 & ev$PEP<=1, ] # >>> 573,970
# Therefore we are separating the results by their protein label.

ev$protLabel <- substr(ev$proteinId1, 1, 3)
ev.con <- subset(ev, ev$protLabel=="CON")
ev.rev <- subset(ev, ev$protLabel=="REV")
ev.psm <- subset(ev, ev$protLabel=="sp|" & ev$PEP<=1) # Only choosing PEP<=1

# Finding the correct sequences from the fasta file:
fasta <- paste(readLines("C:/Users/nslavov/Dropbox/Genetics/Prot Databases/human.fasta"), collapse="")
fasta <- gsub(">", "\n>", fasta) # Each protein in a separate line
writeLines(fasta, "C:/Users/nslavov/Desktop/Percolator/FASTA.txt")

#install.packages("stringr")
library(stringr) # For string matching
#install.packages("Kmisc")
library(Kmisc) # For Reverse

# FOR PSMs:
ev.psm$protLabel <- NULL
ev.psm$proteinId1 <- gsub("[|]", "[|]", ev.psm$proteinId1)
ev.psm$Modified.sequence <- gsub("_", "", ev.psm$Modified.sequence)
protlist <- unique(ev.psm$proteinId1)
ev.psm.seq <- data.frame(id=character(), label=numeric(), ScanNr=numeric(), RT=numeric(), dM=numeric(), PEP=numeric(), PepLen=numeric(), Mass=numeric(),             
                         m.z=numeric(), Score=numeric(), dScore=numeric(), Resolution=numeric(), enzInt=numeric(), n.matches=numeric(), peptide=numeric(),
                         proteinId1=character(), Modified.sequence=character(), Charge1=numeric(), Charge2=numeric(), Charge3=numeric(), Charge4=numeric(),
                         Charge5=numeric(), Charge6=numeric(), Seq=character())
ev.psm$Seq <- NA
pb <- winProgressBar(title="Fixing PSM sequences ...", label="0% done", min=0, max=100, initial=0) # Progress Bar
c <- 0
for (i in protlist) {
  c <- c+1 # Progress Bar
  protSeq <- str_match(fasta, paste(">", i, ".*", sep=""))  # ">sp[|]P04637[|]P53_HUMAN.*"
  ev.psm.p <- subset(ev.psm, ev.psm$proteinId1==i)
  row.names(ev.psm.p) <- NULL
  
  for (ii in 1:nrow(ev.psm.p)) {
    pep <- str_match(protSeq, paste("[A-Z]", ev.psm.p$peptide[ii], "[A-Z]", sep="")) # "[A-Z]ACDGVVHTPAEPTGDSR[A-Z]" 
    if (is.na(pep)){
      pep <- str_match(protSeq, paste(ev.psm.p$peptide[ii], "[A-Z]", sep=""))
      if (is.na(pep)){
        pep <- str_match(protSeq, paste("[A-Z]", ev.psm.p$peptide[ii], sep=""))
        ev.psm.p$Seq[ii] <- paste(substr(pep, 1,1), ".", ev.psm.p$Modified.sequence[ii], ".-", sep="")
      }
      ev.psm.p$Seq[ii] <- paste("-.", ev.psm.p$Modified.sequence[ii], ".", substr(pep, nchar(pep), nchar(pep)), sep="")
    }
    ev.psm.p$Seq[ii] <- paste(substr(pep, 1,1), ".", ev.psm.p$Modified.sequence[ii], ".", substr(pep, nchar(pep), nchar(pep)), sep="")
  }
  ev.psm.seq <- rbind(ev.psm.seq, ev.psm.p)
  
  info <- sprintf("%d%% done", round((c/length(protlist))*100)) # Progress Bar
  setWinProgressBar(pb, c/(length(protlist))*100, label=info) # Progress Bar
}
close(pb) # Progress Bar

write.table(ev.psm.seq, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.psm.seq-v1.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)
ev.psm.seq <- read.delim('C:/Users/nslavov/Desktop/Percolator/MQ files/ev.psm.seq-v1.txt', header = TRUE, stringsAsFactors = FALSE)

#Removing the peptides which were not found in fasta file:
ev.psm.seq$Seq2 <- gsub("NA\\.", NA, ev.psm.seq$Seq)
ev.psm.seq$Seq <- ev.psm.seq$Seq2
ev.psm.seq.not.found <- ev.psm.seq[is.na(ev.psm.seq$Seq2), ]
ev.psm.seq <- ev.psm.seq[!is.na(ev.psm.seq$Seq2), ]
ev.psm.seq$Seq2 <- NULL

ev.psm.seq$peptide <- gsub("\\([a-z]{2})", "*", ev.psm.seq$Seq)
ev.psm.seq$proteinId1 <- gsub("\\[.]", "|", ev.psm.seq$proteinId1)
write.table(ev.psm.seq, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.psm.seq-all.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

ev.psm.seq$Modified.sequence <- NULL
ev.psm.seq$Seq <- NULL
ev.psm.seq <- subset(ev.psm.seq, select = c(1:12, 17:22, 13:16))
write.table(ev.psm.seq, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.psm.perc.ready.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

#FOR REVERSE
ev.rev$protLabel <- NULL
ev.rev$Modified.sequence <- gsub("_", "", ev.rev$Modified.sequence)
res.i <- c(unique(substr(ev.psm.seq$peptide, 1,1)), "-")

for (i in 1:nrow(ev.rev)) { # randomly assigning residues before/after rev sequences
  ev.rev$Seq[i] <- paste(sample(res.i, 1), ".", ev.rev$Modified.sequence[i], ".", sample(res.i, 1), sep="")
}
ev.rev$peptide <- gsub("\\([a-z]{2})", "*", ev.rev$Seq)
write.table(ev.rev, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.rev-all.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)
ev.rev$Modified.sequence <- NULL
ev.rev$Seq <- NULL
ev.rev <- subset(ev.rev, select = c(1:12, 17:22, 13:16))
write.table(ev.rev, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.rev.perc.ready.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

#APPENDING REVERSE AND PSM SEQUENCES:
ev.perc <- rbind(ev.rev, ev.psm.seq)
ev.perc <- ev.perc[sample(nrow(ev.perc)), ] # shuffling the rows
row.names(ev.perc) <- NULL
write.table(ev.perc, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.ready-wRT.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)


#+------------------------------+--------------------------------------------------+
#| PREPARING THE EVIDENCE FILE: | enzN/enzC/ln(numProt)/ln(numPep)/ln(pepSite)/dRT |
#+------------------------------+--------------------------------------------------+

ev.perc.dRT <- read.delim('C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.ready-wRT.txt', header = TRUE, stringsAsFactors = FALSE)
shift.coeffs <- read.delim('C:/Users/nslavov/Desktop/Retention Time Library/shift.coeffs-MQ-Com(0.999-0.9).txt', header = TRUE, stringsAsFactors = FALSE)
rt.lib <- read.delim('C:/Users/nslavov/Desktop/Retention Time Library/RT.library.w-Mean & Median & SD-MQ-Com(0.999-0.9).txt', header = TRUE, stringsAsFactors = FALSE)
ev <- read.delim('C:/Users/nslavov/Desktop/Retention Time Library/evidence.txt', header = TRUE, stringsAsFactors = FALSE)
ev <- ev[ev$Raw.file %in% shift.coeffs[!is.na(shift.coeffs$coeff), "experiment"],
         c('Raw.file', 'Sequence', 'PEP', 'Leading.razor.protein')]

ev.perc.dRT$raw <- gsub("_ID.*", "", ev.perc.dRT$id)
ev.perc.dRT <- ev.perc.dRT [ev.perc.dRT $raw %in% shift.coeffs[!is.na(shift.coeffs$coeff), "experiment"], ] # keeping only the experiments w coeffs
ev.perc.dRT$Seq <- substr(ev.perc.dRT$peptide, 3, nchar(ev.perc.dRT$peptide)-2)
ev.perc.dRT$Seq <- gsub("\\*", "", ev.perc.dRT$Seq)



# Calculating enzN / enzC:
ev.perc.dRT$enzN <- substr(ev.perc.dRT$peptide, 1, 1)
ev.perc.dRT$enzC <- substr(ev.perc.dRT$peptide, nchar(ev.perc.dRT$peptide)-2, nchar(ev.perc.dRT$peptide)-2)
ev.perc.dRT[(ev.perc.dRT$enzN=="K" | ev.perc.dRT$enzN=="R"), "enzN"] <- 1
ev.perc.dRT[!ev.perc.dRT$enzN==1, "enzN"] <- 0
ev.perc.dRT[(ev.perc.dRT$enzC=="K" | ev.perc.dRT$enzC=="R"), "enzC"] <- 1
ev.perc.dRT[!ev.perc.dRT$enzC==1, "enzC"] <- 0
#Be careful some of the residues here are *! but because I know I don't have a modification after
#R or K, I didn't go further to find the residue before *!
ev.perc.dRT <- subset(ev.perc.dRT, select = c(1:18, 25:26, 19:24))



# Calculating ln(numProt) / ln(numPep) / ln(pepSite):
ev.perc.dRT.1 <- data.frame(id=character(), label=numeric(), ScanNr=numeric(), RT=numeric(), dM=numeric(), PEP=numeric(),
PepLen=numeric(), Mass=numeric(), m.z=numeric(), Score=numeric(), dScore=numeric(), Resolution=numeric(), Charge1=numeric(),
Charge2=numeric(), Charge3=numeric(), Charge4=numeric(), Charge5=numeric(), Charge6=numeric(), enzN=numeric(), enzC=numeric(),
enzInt=numeric(), n.matches=numeric(), peptide=character(), proteinId1=character(), raw=character(), Seq=character())

exps <- unique(ev$Raw.file)

for (i in exps) {
  ev.count <- subset(ev, ev$Raw.file==i)
  ev.perc.exp <- subset(ev.perc.dRT, ev.perc.dRT$raw==i)
  
  #ln(numProt):
  prot.count <- as.data.frame(table(ev.count$Leading.razor.protein))
  prot.count$Freq <- log(prot.count$Freq, base=exp(1))
  ev.perc.exp$lnNumProt <- prot.count$Freq[match(ev.perc.exp$proteinId1, prot.count$Var1)]
  
  #ln(numPep):
  pep.count <- as.data.frame(table(ev.count$Sequence))
  pep.count$Freq <- log(pep.count$Freq, base=exp(1))
  ev.perc.exp$lnNumPep <- pep.count$Freq[match(ev.perc.exp$Seq, pep.count$Var1)]
  
  #ln(pepSite):
  prots <- unique(ev.count$Leading.razor.protein)
  proPep <- data.frame(protein=as.character(rep(NA, length(prots))), lnPepSite=0, stringsAsFactors=F)
  c <- 0
  for (ii in prots) {
    c <- c+1
    pep.site <- subset(ev.count, ev.count$Leading.razor.protein==ii)
    proPep$protein[c] <- ii
    proPep$lnPepSite[c] <- log(length(unique(pep.site$Sequence)), base=exp(1))
  }
  ev.perc.exp$lnPepSite <- proPep$lnPepSite[match(ev.perc.exp$proteinId1, proPep$protein)]
  
  #Appending experiments:
  ev.perc.dRT.1 <- rbind(ev.perc.dRT.1, ev.perc.exp)
}
ev.perc.dRT <- ev.perc.dRT.1
ev.perc.dRT <- subset(ev.perc.dRT, select = c(1:18, 27:29, 19:26))



# Separating REV and PSM and calculating dRT:
ev.perc.dRT.p <- subset(ev.perc.dRT, ev.perc.dRT$label==1)

# We only want shared peptides with the library:
ev.perc.dRT.p <- ev.perc.dRT.p[ev.perc.dRT.p$Seq %in% rt.lib$peptides, ]

# PSM:
ev.perc.dRT.p$rt.lib <- rt.lib$rt.median[match(ev.perc.dRT.p$Seq, rt.lib$peptides)]
ev.perc.dRT.p$intercept <- shift.coeffs$intercept[match(ev.perc.dRT.p$raw, shift.coeffs$experiment)]
ev.perc.dRT.p$coeff <- shift.coeffs$coeff[match(ev.perc.dRT.p$raw, shift.coeffs$experiment)]
ev.perc.dRT.p$RT.shift <- ev.perc.dRT.p$intercept + (ev.perc.dRT.p$coeff * ev.perc.dRT.p$RT)
ev.perc.dRT.p$dRT <- abs(ev.perc.dRT.p$RT.shift-ev.perc.dRT.p$rt.lib)

#write.table(ev.perc.dRT.p, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.dRT.p-MQ-Com(0.999-0.9).txt', sep = '\t',
#            row.names = FALSE, quote = FALSE)

ev.perc.dRT.p <- ev.perc.dRT.p[,-(28:33)]

# REV:
# Assigning rt.lib to decoys and calculating dRT for them:
ev.perc.dRT.r <- ev.perc.dRT[ev.perc.dRT$label==-1, ]
ev.perc.dRT.r$intercept <- shift.coeffs$intercept[match(ev.perc.dRT.r$raw, shift.coeffs$experiment)]
ev.perc.dRT.r$coeff <- shift.coeffs$coeff[match(ev.perc.dRT.r$raw, shift.coeffs$experiment)]
ev.perc.dRT.r$RT.shift <- ev.perc.dRT.r$intercept + (ev.perc.dRT.r$coeff * ev.perc.dRT.r$RT)
ev.perc.dRT.r$rt.lib <- sample(rt.lib$rt.median, nrow(ev.perc.dRT.r), replace=T)
ev.perc.dRT.r$dRT <- abs(ev.perc.dRT.r$RT.shif - ev.perc.dRT.r$rt.lib)

write.table(ev.perc.dRT.r, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.dRT.r-MQ-Com(0.999-0.9).txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

ev.perc.dRT.r <- ev.perc.dRT.r[,-(28:33)]

# Combining decoys and PSMs and preparing PIN file:
ev.perc.dRTed <- rbind(ev.perc.dRT.r, ev.perc.dRT.p)
row.names(ev.perc.dRTed) <- NULL
ev.perc.dRTed <- ev.perc.dRTed[sample(nrow(ev.perc.dRTed)), ]
row.names(ev.perc.dRTed) <- NULL

write.table(ev.perc.dRTed, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.dRT+RT(FINAL)-MQ-Com(0.999-0.9)-enzNC-lnPepProSite.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)


#+-----------------------------------------------+------------------------------+
#| PREPARING PERC INPUT FOR DIFFERENT SCENARIOS: | wo RT/RT/dRT/dM*dRT/exp(dRT) |
#+-----------------------------------------------+------------------------------+
ev.perc.dRTed <- read.delim('C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.dRT+RT(FINAL)-MQ-Com(0.999-0.9).txt',
                            header = TRUE, stringsAsFactors = FALSE)
ev.perc.dRTed <- read.delim('C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.dRT+RT(FINAL)-MQ-Com(0.999-0.9)-enzNC.txt',
                            header = TRUE, stringsAsFactors = FALSE)
ev.perc.dRTed <- read.delim('C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.dRT+RT(FINAL)-MQ-Com(0.999-0.9)-enzNC-lnPepProSite.txt',
                            header = TRUE, stringsAsFactors = FALSE)
# +----------+
# | w/o dRT: |
# +----------+
A0 <- ev.perc.dRTed[grep("160902A_NC_set#30A_180min_200NL_60C_N3_250mLC1B", ev.perc.dRTed$id), ]
A0$RT <- NULL
A0$dRT <- NULL
A0$Resolution <- NULL
write.table(A0, 'C:/Users/nslavov/Desktop/Percolator/mq-30A-wo-RT.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)

A0.per <- read.delim('C:/Users/nslavov/Desktop/Percolator/mq-30A-mydRT.txt', header = TRUE, stringsAsFactors = FALSE)
A0.per2 <- A0.per[A0.per$posterior_error_prob<0.02,]

for (i in 1:nrow(A0.per)) {
  A0.per$ev.ID[i] <- str_match(A0.per$PSMId[i], "ID(\\d+)")[2]
}

for (i in 1:nrow(A0)) {
  A0$ev.ID[i] <- str_match(A0$id[i], "ID(\\d+)")[2]
}

rev <- subset(A0, A0$label==-1)
forw <- subset(A0, A0$label==1 & A0$PEP<0.02)

den.for <- density(forw$dRT)
den.rev <- density(rev$dRT)

A0.per$dRT <- A0$dRT[match(A0.per$ev.ID, A0$ev.ID)]
A0.per$PEP <- A0$PEP[match(A0.per$ev.ID, A0$ev.ID)]

A0.per$Tr <- approx(den.for$x, den.for$y, xout=A0.per$dRT)$y
A0.per$Fa <- approx(den.rev$x, den.rev$y, xout=A0.per$dRT)$y
A0.per$PEP.new <- (A0.per$posterior_error_prob*A0.per$Fa)/((A0.per$posterior_error_prob*A0.per$Fa)+((1-A0.per$posterior_error_prob)*A0.per$Tr))

A0.per$Seq <- substr(A0.per$peptide, 3, nchar(A0.per$peptide)-2)
A0.per$Seq <- gsub("\\*", "", A0.per$Seq)

A0.per$Seq <- substr(A0.per$peptide, 3, nchar(A0.per$peptide)-2)
A0.per$Seq <- gsub("\\*", "", A0.per$Seq)

## Preparing all experiments:
exps <- unique(stats.exp$experiment)

for (i in exps) {
  path <- paste('C:/Users/nslavov/Desktop/Percolator/all_exps/all.features.wo.RT/input/', i, '.pin', sep="")
  A0 <- ev.perc.dRTed[grep(i, ev.perc.dRTed$id), ]
  A0$RT <- NULL
  A0$dRT <- NULL
  A0$Resolution <- NULL
  write.table(A0, path, sep = '\t', row.names = FALSE, quote = FALSE)
}

## Running Percolator on all experiments:
exps <- unique(stats.exp$experiment)

for (i in exps) {
  commands <- paste("percolator C:\\Users\\nslavov\\Desktop\\Percolator\\all_exps\\all.features.wo.RT\\input\\", i,
                    ".pin -r C:\\Users\\nslavov\\Desktop\\Percolator\\all_exps\\all.features.wo.RT\\output\\", i, ".txt", sep="")
  system(commands, intern = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = TRUE,
         minimized = FALSE, invisible = TRUE)
}

## Combining all of the experiments into one file:
all.perc <- data.frame(PSMId=character(), score=numeric(), q.value=numeric(), posterior_error_prob=numeric(),
                       peptide=character(), proteinIds=character(), stringsAsFactors=F)

exp.names <- data.frame(name=list.files('C:/Users/nslavov/Desktop/Percolator/all_exps/all.features.wo.RT/output'), stringsAsFactors=F)


for (i in exp.names$name) {
  path <- paste("C:/Users/nslavov/Desktop/Percolator/all_exps/all.features.wo.RT/output/", i, sep="")
  exp.i <- read.delim(path, header=T, stringsAsFactors=F, sep='\t')
  all.perc <- rbind(all.perc, exp.i)
}

all.perc$peptide <- as.character(all.perc$peptide)
all.perc$Seq <- substr(all.perc$peptide, 3, nchar(all.perc$peptide)-2)
all.perc$Seq <- gsub("\\*", "", all.perc$Seq)

write.table(all.perc, 'C:/Users/nslavov/Desktop/Percolator/all_exps/all.features.wo.RT/output/all.exps-all.feats.wo.RT.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

## Bayesian Updating Perc Results:
for (i in 1:nrow(all.perc)) {
  all.perc$ev.ID[i] <- str_match(all.perc$PSMId[i], "ID(\\d+)")[2]
}
for (i in 1:nrow(ev.perc.dRTed)) {
  ev.perc.dRTed$ev.ID[i] <- str_match(ev.perc.dRTed$id[i], "ID(\\d+)")[2]
}
write.table(ev.perc.dRTed, 'C:/Users/nslavov/Desktop/Percolator/MQ files/ev.perc.dRT+RT(FINAL)-MQ-Com(0.999-0.9)-enzNC-lnPepProSite+ev.ID extracted.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

all.perc.ed <- data.frame(PSMId=character(), score=numeric(), q.value=numeric(), posterior_error_prob=numeric(),
                          peptide=character(), proteinIds=character(), Seq=character(), ev.ID=numeric(), dRT=numeric(),
                          PEP=numeric(), Tr=numeric(), Fa=numeric(), PEP.new=numeric())
for (i in exps) {
  A <- all.perc[grep(i, all.perc$PSMId), ]
  A0 <- ev.perc.dRTed[grep(i, ev.perc.dRTed$id), ]
  A$dRT <- A0$dRT[match(A$ev.ID, A0$ev.ID)]
  
  if (length(A$PSMId)==0){ # Because all of the exps are not processed by Percolator (95/118 are processed).
    next
  }
  
  rev <- subset(A0, A0$label==-1)
  #forw <- subset(A0, A0$label==1 & A0$PEP<0.02) Defining correct hits based on MQ or Percolator
  forw <- subset(A, A$posterior_error_prob<0.02) # Based on Percolator
  
  
  den.for <- density(forw$dRT)
  den.rev <- density(rev$dRT)
  
  #A$dRT <- A0$dRT[match(A$ev.ID, A0$ev.ID)]
  A$PEP <- A0$PEP[match(A$ev.ID, A0$ev.ID)]
  
  A$Tr <- approx(den.for$x, den.for$y, xout=A$dRT)$y
  A$Fa <- approx(den.rev$x, den.rev$y, xout=A$dRT)$y
  A$PEP.new <- (A$posterior_error_prob*A$Fa)/((A$posterior_error_prob*A$Fa)+((1-A$posterior_error_prob)*A$Tr))

  all.perc.ed <- rbind(A, all.perc.ed)  
}
write.table(all.perc.ed, 'C:/Users/nslavov/Desktop/Percolator/all_exps/all.features.wo.RT/output/all.exps-all.feats.wo.RT+Bayesian.dRT.Update+Perc.forw.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

# +-------------------------------+
# | Using K option of Percolator: |
# +-------------------------------+
A1 <- ev.perc.dRTed[grep("160902A_NC_set#30A_180min_200NL_60C_N3_250mLC1B", ev.perc.dRTed$id), ]
A1$dRT <- NULL
A1$Resolution <- NULL
write.table(A1, 'C:/Users/nslavov/Desktop/Percolator/mq-h50-percRT.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)

A02 <- subset(A, A$PEP<=0.02)
A02$Seq <- substr(A02$peptide, 3, nchar(A02$peptide)-2)
A02$Seq <- gsub("\\*", "", A02$Seq)
pep <- unique(A02$Seq)
FDR <- mean(A02$PEP)


# +---------------+
# | Using my dRT: |
# +---------------+
A2 <- ev.perc.dRTed[grep("160808A_NC_set19A_180min_200NL_targeted", ev.perc.dRTed$id), ]
A2$RT <- A2$dRT
A2$dRT <- NULL
A2$Resolution <- NULL
colnames(A2)[4] <- "dRT"
write.table(A2, 'C:/Users/nslavov/Desktop/Percolator/mq-19A-mydRT.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)

# Creating input for all of the experiments:

exps <- unique(stats.exp$experiment)

for (i in exps) {
  path <- paste('C:/Users/nslavov/Desktop/Percolator/all_exps/', i, '.pin', sep="")
  A2 <- ev.perc.dRTed[grep(i, ev.perc.dRTed$id), ]
  A2$RT <- A2$dRT
  A2$dRT <- NULL
  A2$Resolution <- NULL
  colnames(A2)[4] <- "dRT"
  write.table(A2, path, sep = '\t', row.names = FALSE, quote = FALSE)
}

A2.per <- read.delim('C:/Users/nslavov/Desktop/Percolator/mq-19A-mydRT.txt', header = TRUE, stringsAsFactors = FALSE)
A2.per$Seq <- substr(A2.per$peptide, 3, nchar(A2.per$peptide)-2)
A2.per$Seq <- gsub("\\*", "", A2.per$Seq)
A2.per2 <- A2.per[A2.per$posterior_error_prob<0.02,]
FDR <- mean(A2.per2$posterior_error_prob) 


# +------------------------+
# | Using my dRT & dRT*dM: |
# +------------------------+
A2 <- ev.perc.dRTed[grep("160801A_NC_set19B_180min_50ID_100NL_from100ul", ev.perc.dRTed$id), ]
A2$RT <- A2$dRT
A2$dRT <- NULL
A2$Resolution <- NULL
colnames(A2)[4] <- "dRT"
A2$dRdM <- A2$dRT * A2$dM
A2 <- subset(A2, select = c(1:5, 22, 6:21))

write.table(A2, 'C:/Users/nslavov/Desktop/Percolator/mq-h50-mydRTdM.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)


# +-------------+
# | w exp(dRT): |
# +-------------+
A3 <- ev.perc.dRTed[grep("160801A_NC_set19B_180min_50ID_100NL_from100ul", ev.perc.dRTed$id), ]
A3$RT <- exp(A3$dRT)
A3$dRT <- NULL
A3$Resolution <- NULL
colnames(A3)[4] <- "dRT"
write.table(A3, 'C:/Users/nslavov/Desktop/Percolator/mq-h50-exp(mydRT).pin', sep = '\t',
            row.names = FALSE, quote = FALSE)

A3.per <- read.delim('C:/Users/nslavov/Desktop/Percolator/mq-h50-mydRT-enzNC.txt', header = TRUE, stringsAsFactors = FALSE)
A3.per2 <- A2.per[A3.per$q.value<0.02,]


# +------------+
# | w ln(dRT): |
# +------------+
A4 <- ev.perc.dRTed[grep("160801A_NC_set19B_180min_50ID_100NL_from100ul", ev.perc.dRTed$id), ]
A4$RT <- log(A4$dRT, base = exp(1))
A4$dRT <- NULL
A4$Resolution <- NULL
colnames(A4)[4] <- "dRT"
write.table(A4, 'C:/Users/nslavov/Desktop/Percolator/mq-h50-ln(mydRT).pin', sep = '\t',
            row.names = FALSE, quote = FALSE)

A3.per <- read.delim('C:/Users/nslavov/Desktop/Percolator/mq-h50-mydRT-enzNC.txt', header = TRUE, stringsAsFactors = FALSE)
A3.per2 <- A2.per[A3.per$q.value<0.02,]
FDR <- mean(A3.per2$posterior_error_prob)


# +------------------------------------+
# | Preparing the whole evidence file: |
# | my dRT:                            |
# +------------------------------------+
A <- ev.perc.dRTed
A$RT <- A$dRT
A$dRT <- NULL
A$Resolution <- NULL
colnames(A)[4] <- "dRT"
A$ScanNr <- 1
write.table(A, 'C:/Users/nslavov/Desktop/Percolator/evidence-mydRT.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)

A.per <- read.delim('C:/Users/nslavov/Desktop/Percolator/evidence-mydRT-enzNC.txt', header = TRUE, stringsAsFactors = FALSE)
A.per2 <- A.per[A.per$posterior_error_prob<0.02,]


# +---------------------------------------------------------+
# | Using Max Quant features for PEP calculation + my dRT   |
# | Andromeda Score/Peptide Length/Charge State             |
# | No. missed cleavages/# of variable modifications        |
# | We do the searches w/o var mods                         |
# | We need the evidence file for # of missed cleavages     |
# +---------------------------------------------------------+
library(stringr)
ev <- read.delim('C:/Users/nslavov/Desktop/Retention Time Library/evidence.txt', header = TRUE, stringsAsFactors = FALSE)

A2 <- ev.perc.dRTed[grep("160902A_NC_set#30A_180min_200NL_60C_N3_250mLC1B", ev.perc.dRTed$id), ]
A2$RT <- A2$dRT
A2$dRT <- NULL
A2$Resolution <- NULL
A2$dM <- NULL
A2$PEP <- NULL
A2$Mass <- NULL
A2$m.z <- NULL
A2$dScore <- NULL
A2$lnNumProt <- NULL
A2$lnNumPep <- NULL
A2$lnPepSite <- NULL
A2$enzN <- NULL
A2$enzC <- NULL
A2$enzInt <- NULL
A2$enzInt <- NULL
A2$n.matches <- NULL

for (i in 1:nrow(A2)) {
  A2$ev.ID[i] <- str_match(A2$id[i], "ID(\\d+)")[2]
}

A2$n.missed.cleavages <- ev$Missed.cleavages[match(A2$ev.ID, ev$id)]
A2 <- subset(A2, select = c(1:12, 16, 13:14))
colnames(A2)[4] <- "dRT"
write.table(A2, 'C:/Users/nslavov/Desktop/Percolator/mq-30A-onlyMQ+mydRT.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)
##
A2.per <- read.delim('C:/Users/nslavov/Desktop/Percolator/mq-30A-onlyMQ+mydRT.txt', header = TRUE, stringsAsFactors = FALSE)
A2.per2 <- A2.per[A2.per$posterior_error_prob<0.02,]























## Other test codes:
A <- ev.t[grep("160801A_NC_set19B_180min_50ID_100NL_from100ul", ev.t$id), ]
A <- A[!is.na(A$dM),]
A <- A[!is.na(A$Resolution),]
write.table(A, 'C:/Users/nslavov/Desktop/Percolator/mq-h50-perc.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)
A2 <- A[A$PEP<0.05, ]
FDRmq <- mean(A2$PEP)


A.per <- read.delim('C:/Users/nslavov/Desktop/Percolator/mq-h50-perc.txt', header = TRUE, stringsAsFactors = FALSE)
A.per2 <- A.per[A.per$q.value<0.01,]

h50.pin$peptide <- comet.pin$peptide[match(h50.pin$ScanNr, comet.pin$ScanNr)]
h50.pin <- h50.pin[,-(28:31)]
h50.pin.RT <- h50.pin[, c("id", "label", "ScanNr", "dRT", "lnrSp", "deltLCn", "deltCn", "lnExpect", "Xcorr", "Sp", "IonFrac", "Mass", "PepLen",
                          "Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "Charge6", "enzN", "enzC", "enzInt", "lnNumSP", "dM", "absdM",
                          "peptide", "proteinId1")] #With my dRTs



ev02 <- ev.perc.dRT[ev.perc.dRT$PEP<=0.05, ]
Z <- length(unique(ev02$Seq))

# Calculating qValue:
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library("qvalue")

ev <- read.delim('C:/Users/nslavov/Desktop/Retention Time Library/evidence.txt', header = TRUE, stringsAsFactors = FALSE)
ev <- ev[, c('Raw.file', 'Sequence', 'PEP', 'Retention.time', 'Scan.index')]









