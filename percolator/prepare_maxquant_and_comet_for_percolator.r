# A huge chunk of the following code prepares mq and comet individually for feeding into
# percolator. Preparation of Comet includes taking the Reten times from the mzXML raw files.
# or calculating the dRT's by matching the comet's result with the rt.lib.


ms2 <- read.delim('C:/Users/nslavov/Desktop/New MSMS/msms.txt', header = TRUE, stringsAsFactors = FALSE)
comet.pin <- read.delim('C:/Users/nslavov/Desktop/Percolator/h50.pin', header=T)
ms2.exp <- subset(ms2, ms2$Raw.file=="160801A_NC_set19B_180min_50ID_100NL_from100ul")
Z <- comet.pin[comet.pin$ScanNr %in% ms2.exp$Scan.number, ]
Z <- Z[, c("id", "label", "ScanNr", "peptide", "proteinId1")]
ms2.exp <- ms2.exp[ms2.exp$Scan.number %in% Z$ScanNr, ]
ms2.exp <- ms2.exp[, c("Raw.file", "Scan.number", "Sequence", "Reverse", "PEP")]

ms2sc <- read.delim('C:/Users/nslavov/Desktop/New MSMS/msmsScans.txt', header = TRUE, stringsAsFactors = FALSE)
ms2sc <- ms2sc[, c("Raw.file", "Scan.number", "Retention.time")]
ms2sc.exp <- ms2sc[ms2sc$Raw.file=="160801A_NC_set19B_180min_50ID_100NL_from100ul", ]
rownames(ms2sc.exp) <- NULL

M <- unique(ms2sc.exp$Scan.number)


# Parsing the mzXML to get the Scan Information:
source("https://bioconductor.org/biocLite.R")
biocLite("mzR")
library("mzR")

ms <- openMSfile("C:/Users/nslavov/Desktop/h50.mzXML")
Scans <- header(ms) #This file has the information in the mzXML file
Scans$ScanNo <- rownames(Scans)
#-------------------------------------------

install.packages("XML")
install.packages("methods")
library("XML")
library("methods")
result <- xmlParse("C:/Users/nslavov/Desktop/Percolator/output.xml")

xmldataframe <- xmlToDataFrame("C:/Users/nslavov/Desktop/Percolator/output.xml")
print(xmldataframe)

mss <- openMSfile("C:/Users/nslavov/Desktop/Percolator/output.xml")

#-------------------------------------------

#PERCOLATOR: 
perc.tab <- read.delim('C:/Users/nslavov/Desktop/Percolator/output.txt', header=T)

perc.tab02 <- subset(perc.tab, perc.tab$q.value<0.02)
Z <- unique(perc.tab02$proteinIds)
FDR <- mean(perc.tab02$posterior_error_prob)

ms2.exp02 <- subset(ms2.exp, ms2.exp$PEP<.03)
ms2.exp02 <- ms2.exp02[ms2.exp02$Sequence %in% rt.lib$peptides, ]
Zmq <- unique(ms2.exp02$Sequence)
FDRmq <- mean(ms2.exp02$PEP)

#-------------------------------------------
# PREPARING PERCULATOR INPUT:

row.names(ms2.exp) <- NULL

ms2sc.exp.pin <- ms2.exp[,c("Raw.file", "Reverse", "Scan.number", "PEP", "Score", "Delta.score", "Mass.Error..ppm.", "Sequence", "Proteins")]
ms2sc.exp.pin$Raw.file <- paste("h50", rownames(ms2sc.exp.pin), sep = "-")
colnames(ms2sc.exp.pin) <- c("id", "label", "ScanNr", "PEP", "Score", "Delta.score", "Mass.Error.ppm", "peptide", "proteinId1")
ms2sc.exp.pin[!ms2sc.exp.pin$label=="+", "label"] <- 1
ms2sc.exp.pin[ms2sc.exp.pin$label=="+", "label"] <- -1

write.table(ms2sc.exp.pin, 'C:/Users/nslavov/Desktop/Percolator/mq-h50.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)


#-------------------------------------------------------------------------------
# PROCESSING THE COMET RESULTS TO ADD DELTA.RTs
comet.pin <- read.delim('C:/Users/nslavov/Desktop/Percolator/h50.pin', header=T)
rt.lib <- read.delim('C:/Users/nslavov/Desktop/Retention Time Library/RT.library.w-Mean & Median & SD.txt', header = TRUE, stringsAsFactors = FALSE)
shift.coeffs <- read.delim('C:/Users/nslavov/Desktop/Retention Time Library/shift.coeffs.txt', header = TRUE, stringsAsFactors = FALSE)

# Cleaning sequences:
co.ed <- comet.pin
co.ed$peptide <- as.character(co.ed$peptide)
co.ed$peptide <- gsub("\\*", "", co.ed$peptide)
co.ed$peptide <- substr(co.ed$peptide, 3, nchar(co.ed$peptide)-2)

# Getting RTs for the peptides:
library("mzR")
h50 <- openMSfile("C:/Users/nslavov/Desktop/h50.mzXML")
h50 <- header(h50) #This has the information in the mzXML file
co.ed$RT <- h50$retentionTime[match(co.ed$ScanNr, h50$seqNum)]/60

# Shifting the RTs:
co.ed$coeff <- shift.coeffs[shift.coeffs$experiment=="160801A_NC_set19B_180min_50ID_100NL_from100ul", "coeff"]
co.ed$intercept <- shift.coeffs[shift.coeffs$experiment=="160801A_NC_set19B_180min_50ID_100NL_from100ul", "intercept"]
co.ed$RT.shif <- co.ed$intercept + (co.ed$coeff * co.ed$RT)

# Keeping the PSMs in the library and calculating deltaRT for them:
co.ed.psm <- co.ed[co.ed$label==1, ]
co.ed.psm <- co.ed.psm[co.ed.psm$peptide %in% rt.lib$peptides, ]

co.ed.psm$rt.lib <- rt.lib$rt.median[match(co.ed.psm$peptide, rt.lib$peptides)]
co.ed.psm$dRT <- abs(co.ed.psm$RT.shif - co.ed.psm$rt.lib)
row.names(co.ed.psm) <- NULL

# Assigning rt.lib to decoys and calculating dRT for them:
co.ed.dec <- co.ed[co.ed$label==-1, ]
co.ed.dec <- co.ed.dec[sample(nrow(co.ed.dec)),]
row.names(co.ed.dec) <- NULL
del <- seq(1, nrow(co.ed.dec), 2) #removing half of decoys randomly because they are twice the size of PSMs
co.ed.dec <- co.ed.dec[-del, ]
row.names(co.ed.dec) <- NULL

co.ed.dec$rt.lib <- sample(rt.lib$rt.median, nrow(co.ed.dec), replace=T)
co.ed.dec$dRT <- abs(co.ed.dec$RT.shif - co.ed.dec$rt.lib)

# Combining decoys and PSMs and preparing PIN file:
h50.pin <- rbind(co.ed.dec,co.ed.psm)
row.names(h50.pin) <- NULL
h50.pin <- h50.pin[sample(nrow(h50.pin)), ]
row.names(h50.pin) <- NULL

h50.pin$peptide <- comet.pin$peptide[match(h50.pin$ScanNr, comet.pin$ScanNr)]
h50.pin <- h50.pin[,-(28:31)]
h50.pin.RT <- h50.pin[, c("id", "label", "ScanNr", "dRT", "lnrSp", "deltLCn", "deltCn", "lnExpect", "Xcorr", "Sp", "IonFrac", "Mass", "PepLen",
                          "Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "Charge6", "enzN", "enzC", "enzInt", "lnNumSP", "dM", "absdM",
                          "peptide", "proteinId1")] #With my dRTs

# Without dRT or RTs:
h50.pin <- h50.pin.RT[,-4]

# With RT for letting Percolator calculate its dRTs:
h50.perc.RT <- h50.pin.RT
h50.perc.RT$dRT <- h50$retentionTime[match(h50.perc.RT$ScanNr, h50$seqNum)]/60
h50.perc.RT <- subset(h50.perc.RT, select = c(1:4, 24, 5:23, 25:27))
colnames(h50.perc.RT)[4] <- "RT"

# With dM*dRT feature (My dRTs):
h50.pin.RT.M <- h50.pin.RT
h50.pin.RT.M$dmdRT <- h50.pin.RT.M$dRT * h50.pin.RT.M$dM
h50.pin.RT.M <- subset(h50.pin.RT.M, select = c(1:4, 28, 5:27))
h50.test <- subset(h50.pin, select = c(1:6, 8:11, 22:26))

write.table(h50.pin.RT, 'C:/Users/nslavov/Desktop/Percolator/h50-dRT.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)
write.table(h50.pin, 'C:/Users/nslavov/Desktop/Percolator/h50-wodRT.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)
write.table(h50.perc.RT, 'C:/Users/nslavov/Desktop/Percolator/h50-RT-K.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)
write.table(h50.pin.RT.M, 'C:/Users/nslavov/Desktop/Percolator/h50-dRTdM.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)
write.table(h50.test, 'C:/Users/nslavov/Desktop/Percolator/h50-test.pin', sep = '\t',
            row.names = FALSE, quote = FALSE)

perc.pin.RT <- read.delim('C:/Users/nslavov/Desktop/Percolator/h50-dRT.txt', header=T)
perc.pin.woRT <- read.delim('C:/Users/nslavov/Desktop/Percolator/h50-wodRT.txt', header=T)
perc.pin.KRT <- read.delim('C:/Users/nslavov/Desktop/Percolator/h50-RT-K.txt', header=T)


perc.pin.RT.t <- subset(perc.pin.RT, perc.pin.RT$q.value<0.03)

FDR <- mean(perc.pin.RT.t$posterior_error_prob)

hist(co.ed.dec$dRT)


co.ed.decoy <- co.ed[co.ed$label==-1, ]

#-------------------------------------------------------------------------------


ev <- read.delim('C:/Users/nslavov/Desktop/New MSMS/evidence.txt', header = TRUE, stringsAsFactors = FALSE)
ev.exp <- ev[ev$Raw.file=="160801A_NC_set19B_180min_50ID_100NL_from100ul", ]
ev.exp.t <- ev.exp[ev.exp$PEP<=1,]
ev.exp.t <- ev.exp.t[ev.exp.t$Reverse=="+", ]
ev <- ev[ev$PEP<=1, c('Raw.file', 'Sequence', 'PEP', 'Retention.time', 'id')]

ev.exp.t1 <- ev.exp[!ev.exp$Reverse=="+", ]
ev.exp.t2 <- ev.exp[ev.exp$Reverse=="+", ]

plot(density(ev.exp.t1$Retention.length, na.rm=T), col="green")
lines(density(ev.exp.t2$Retention.length, na.rm=T))

#-------------------------------------------------------------------------------

#=======================================================================================

# FIXING MASS CALCULATIONS AND dM:
ev.perc$Seq <- substr(ev.perc$peptide, 3, nchar(ev.perc$peptide)-2)
ev.perc$nTMT <- str_count(ev.perc$Seq,"K")+1
ev.perc$nOX <- str_count(ev.perc$Seq,"ox")

ev.psm.seq <- read.delim('C:/Users/nslavov/Desktop/Percolator/MQ files/ev.psm.seq-all', header = TRUE, stringsAsFactors = FALSE)
# Fixing Mass and dM:
ev.psm.seq$nTMT <- str_count(ev.psm.seq$Modified.sequence,"K")+1 # +1 is for N-terminus
ev.psm.seq$nOX <- str_count(ev.psm.seq$Modified.sequence,"ox")
ev.psm.seq$nPH <- str_count(ev.psm.seq$Modified.sequence,"ph")
ev.psm.seq$nDE <- str_count(ev.psm.seq$Modified.sequence,"de")
ev.psm.seq$charge <- 0
ev.psm.seq[ev.psm.seq$Charge1==1, "charge"] <- 1
ev.psm.seq[ev.psm.seq$Charge2==1, "charge"] <- 2
ev.psm.seq[ev.psm.seq$Charge3==1, "charge"] <- 3
ev.psm.seq[ev.psm.seq$Charge4==1, "charge"] <- 4
ev.psm.seq[ev.psm.seq$Charge5==1, "charge"] <- 5
ev.psm.seq[ev.psm.seq$Charge6==1, "charge"] <- 6

# MQ has not taken into account the mass of TMT labels and Hydrogen when
# it calculated the Mass. We should just add the mass of TMT and H:
# OX(Oxidation): 15.9949146221 / PH(Phosphorylation): 79.9663304084 / DE(Deamidation (NQ)): 0.9840155848


ev.psm.seq$theoMass <- ev.psm.seq$Mass + ev.psm.seq$charge*1.007276 # + 229.162932*ev.psm.seq$nTMT
ev.psm.seq$obsMass <- ev.psm.seq$m.z * ev.psm.seq$charge
ev.psm.seq$dM.new <- abs(ev.psm.seq$theoMass-ev.psm.seq$obsMass)
ev.psm.seq$dM.new <- ev.psm.seq$dM - (ev.psm.seq$charge*1.007276 + 229.162932*ev.psm.seq$nTMT)/ev.psm.seq$charge








160801A_NC_set19B_180min_50ID_100NL_from100ul



















ms2 <- read.delim('C:/Users/nslavov/Desktop/New MSMS/msms.txt', header = TRUE, stringsAsFactors = FALSE)
ms2.exp <- subset(ms2, ms2$Raw.file=="160801A_NC_set19B_180min_50ID_100NL_from100ul")
ms2.exp.t1 <- ms2.exp[ms2.exp$PEP<=0.05, ]
ms2.exp.t2 <- ms2.exp[ms2.exp$PEP>0.05, ]
ms2.exp.t <- ms2.exp.t[ms2.exp.t$Reverse=="+", ]

ms2.exp.t1 <- ms2.exp[!ms2.exp$Reverse=="+", ]
ms2.exp.t2 <-  ms2.exp[ms2.exp$Reverse=="+", ]

ms2.exp <- ms2.exp[, c("Raw.file", "Reverse", "Scan.number", "Sequence", "PEP", )]


plot(density(ms2.exp.t1$Retention.l, na.rm=T), col="green")
lines(density(ms2.exp.t2$Precursor.Intensity, na.rm=T))

ms2 <- read.delim('C:/Users/nslavov/Desktop/New MSMS/evidence+PEP(1)+DC.txt', header = TRUE, stringsAsFactors = FALSE)











# TEST CODE:
A <- regexpr("\n>sp|Q0P651|ABD18_HUMAN.*$", FASt, perl = T)
B <- gregexpr(">sp|Q0P651|ABD18_HUMAN.*([A-Z]GVSKL[A-Z])", FASt)
A <- grep(">sp[|]P04637[|]P53_HUMAN.*", FASt, perl = T, value = T)
FASt <- gsub(">sp[|]P04637[|]P53_HUMAN.*", A, FASt, perl=T)
writeLines(A, "C:/Users/nslavov/Desktop/Percolator/A.txt")
