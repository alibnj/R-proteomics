msms <- read.delim('C:/Users/nslavov/Desktop/R-PSM vs Duty Cycles/msms.txt', header=TRUE, stringsAsFactors=FALSE)
spec <- read.delim('C:/Inetpub/wwwroot/ISB/data/test3/MaxQuant.VS.Spec.v3.of.TPP-XT/h129-V3-BIN4-pep0.9-13.xls', header=TRUE, stringsAsFactors=FALSE)

msms1 <- subset(msms, Raw.file == '1608329A_NC_set#29B_180min_200NL_60C+N=3_250ms')
msms1 <- subset(msms1, PEP < 0.02)
uni.seq <- as.data.frame(unique(msms1$Sequence))
pep <- as.data.frame(unique(msms1$Sequence))

A <- substr(spec$peptide, 12, nchar(spec$peptide)-2) #extract characters from a cell with start and end point
B <- gsub("[[357.26]]", "", A)
B <- gsub("[[357.2]", "", B)
B <- as.data.frame(B)
Bu <- as.data.frame(unique(B))
colnames(Bu) <- c('A')

AA <- as.data.frame(intersect(Bu$A, uni.seq$`unique(msms1$Sequence)`))

peppro <- read.delim('C:/Inetpub/wwwroot/ISB/data/test3/MaxQuant.VS.Spec.v3.of.TPP-XT/h10-V3-BIN4-pep0.9-3131.xls', header=TRUE, stringsAsFactors=FALSE)
meann <- mean(peppro$probability)
FDR <- 1-meann

hist(peppro$probability)

ggplot(peppro, aes(probability, fill = 'b', colour = 'b')) +
  geom_density(alpha=0.9)+xlab('Peptide Prophet Probability')+
  ylab('Probability')+#labs(title=experiments[e,1])+
  theme(text = element_text(size=20))

uni.seq$is <- Bu$A[match(uni.seq$`unique(msms1$Sequence)`, Bu$A)]

#-------------------------------
#---COMPARING MQ and SpectraST--
#-------------------------------
MSc <- read.csv('C:/Users/nslavov/Desktop/R-LogReg/msmsScan-With PEP.csv', header=TRUE, stringsAsFactors=FALSE)

for (i in 1:14)
{
  mq <- paste("C:/Users/nslavov/Desktop/MQ.vs.Spec/", i, ".txt", sep = "")
  sp <- paste("C:/Users/nslavov/Desktop/MQ.vs.Spec/Spec/ed/", i, ".ipro.098.xls", sep = "")
  ms <- paste("C:/Users/nslavov/Desktop/MQ.vs.Spec/Spec/ed/", i, ".m.s.csv", sep = "")
  m <- read.delim(mq, header=TRUE)
  s <- read.delim(sp, header=TRUE, sep = '\t')
  
  m$Sequence <- as.character(m$Sequence)
  s <- as.data.frame(unique(s$peptide))
  colnames(s) <- c('spec.pep')
  s$spec.pep <- as.character(s$spec.pep)
  
  s$spec.pep <- substr(s$spec.pep, 3, nchar(s$spec.pep)-2) #extract characters from a cell with start and end point
  shared <- as.data.frame(intersect(m$Sequence, s$spec.pep))
  write.csv(shared, ms)
}



#-----------------------------------------------
#---COMPARING MQ and SPECTRAST (V3) OVERLAPS----
#-----------------------------------------------
msms <- read.delim('C:/Users/nslavov/Desktop/R-LogReg/msms.txt', header=TRUE, stringsAsFactors=FALSE)
msms <- subset(msms,select=c('Raw.file', 'Sequence', 'Charge', 'PEP'))

msms[is.na(msms$PEP),] <- NULL
msms2 <- msms[!msms$PEP>1,]
msms2 <- msms2[!msms2$Charge==0,]
msms <- msms2

spec <- read.delim('C:/Inetpub/wwwroot/ISB/data/test3/MaxQuant.VS.Spec.v3.of.TPP-XT/h129-V3-BIN4-pep0.9-13.xls', header=TRUE, stringsAsFactors=FALSE)

msms1 <- subset(msms, Raw.file == '1608329A_NC_set#29B_180min_200NL_60C+N=3_250ms')
msms1 <- subset(msms1, PEP < 0.02)
uni.seq <- as.data.frame(unique(msms1$Sequence))

A <- substr(spec$peptide, 12, nchar(spec$peptide)-2) #extract characters from a cell with start and end point
B <- gsub("[[357.26]]", "", A)
B <- gsub("[[357.2]", "", B)
B <- as.data.frame(B)
Bu <- as.data.frame(unique(B))
colnames(Bu) <- c('A')

AA <- as.data.frame(intersect(Bu$A, uni.seq$`unique(msms1$Sequence)`))
AA.rem <- intersect(Bu$A, uni.seq$`unique(msms1$Sequence)`)

distinct.spec <- as.data.frame(Bu[ ! Bu$A %in% AA.rem, ])
colnames(distinct.spec) <- c('A')

msms3 <- subset(msms, Raw.file == '1608329A_NC_set#29B_180min_200NL_60C+N=3_250ms')
msms3 <- subset(msms3, PEP > 0.02)
uni.seq.3 <- as.data.frame(unique(msms3$Sequence))
colnames(uni.seq.3) <- c('A')
BB <- as.data.frame(intersect(distinct.spec$A, uni.seq.3$A))

#-------------------------------
#---COMET FDR-------------------
#-------------------------------

### DATA CLEANING FOR GETTING THE OVERLAP WITH MAX QUANT ####

comet <- read.delim('C:/Inetpub/wwwroot/ISB/data/comet-cluster/h1.V1-095.pep-095.xls', header=TRUE, stringsAsFactors=FALSE)

msms1 <- subset(msms, Raw.file == '1608329A_NC_set#29B_180min_200NL_60C+N=3_250ms')
msms1 <- subset(msms1, PEP < 0.02)
uni.seq <- as.data.frame(unique(msms1$Sequence))
pep <- as.data.frame(unique(msms1$Sequence))

A <- substr(comet$peptide, 12, nchar(comet$peptide)-2) #extract characters from a cell with start and end point
B <- gsub("[[357.26]]", "", A)
B <- gsub("[[357.2]", "", B)
B <- as.data.frame(B)
Bu <- as.data.frame(unique(B))
colnames(Bu) <- c('A')

AA <- as.data.frame(intersect(Bu$A, uni.seq$`unique(msms1$Sequence)`))

#### FDR ####

peppro <- read.delim('C:/Inetpub/wwwroot/ISB/data/spectrast.vs.comet-V1.library/h129.V1-095.pep-095.xls', header=TRUE, stringsAsFactors=FALSE)
meann <- mean(peppro$probability)
FDR <- 1-meann

#############
 
hist(peppro$probability)

ggplot(peppro, aes(probability, fill = 'b', colour = 'b')) +
  geom_density(alpha=0.9)+xlab('Peptide Prophet Probability')+
  ylab('Probability')+#labs(title=experiments[e,1])+
  theme(text = element_text(size=20))

uni.seq$is <- Bu$A[match(uni.seq$`unique(msms1$Sequence)`, Bu$A)]

#--------------------------------------------------------
#---COMET VS. SPECTRAST (LIB COMET V1)-------------------
#--------------------------------------------------------


#### COMET:
com <- read.delim('C:/Inetpub/wwwroot/ISB/data/spectrast.vs.comet-V1.library/comet/h129.pep.xls', header=TRUE, stringsAsFactors=FALSE)

#---Cleaning:
A <- substr(com$peptide, 12, nchar(com$peptide)-2) #extract characters from a cell with start and end point
B <- gsub("[[357.26]]", "", A)
B <- gsub("[[357.2]", "", B)
B <- gsub("1404]", "", B)
B <- gsub("1600", "", B)
B <- as.data.frame(B)
com.seqs <- as.data.frame(unique(B))
colnames(com.seqs) <- c('A')

#### SpecraST:
spec <- read.delim('C:/Inetpub/wwwroot/ISB/data/spectrast.vs.comet-V1.library/h129.V1-095.pep-095.xls', header=TRUE, stringsAsFactors=FALSE)

#---Cleaning:
A <- substr(spec$peptide, 12, nchar(spec$peptide)-2) #extract characters from a cell with start and end point
B <- gsub("[[357.26]]", "", A)
B <- gsub("[[357.2]", "", B)
B <- gsub("1404]", "", B)
B <- gsub("1600", "", B)
B <- as.data.frame(B)
spec.seqs <- as.data.frame(unique(B))
colnames(spec.seqs) <- c('A')

#### Overlap:
overlap <- as.data.frame(intersect(spec.seqs$A, com.seqs$A))
overlap.rem <- intersect(spec.seqs$A, com.seqs$A)

distinct.spec <- as.data.frame(spec.seqs[ ! spec.seqs$A %in% overlap.rem, ]) #Hits that are only in SpectraST results
colnames(distinct.spec) <- c('A')

#### COMET with P.min=0:
com.0 <- read.delim('C:/Inetpub/wwwroot/ISB/data/spectrast.vs.comet-V1.library/comet.pep0/h129.pep.xls', header=TRUE, stringsAsFactors=FALSE)
com.0 <- com.0[ ! com.0$probability=='[unavailable]' , ]
com.95 <- subset(com.0, probability < 0.95 & probability >= 0.5)
com.50 <- subset(com.0, probability < 0.5)

# Comparison with P<0.95:
A <- substr(com.95$peptide, 12, nchar(com.95$peptide)-2) #extract characters from a cell with start and end point
B <- gsub("[[357.26]]", "", A)
B <- gsub("[[357.2]", "", B)
B <- gsub("1404]", "", B)
B <- gsub("1600", "", B)
B <- as.data.frame(B)
com.95.seqs <- as.data.frame(unique(B))
colnames(com.95.seqs) <- c('A')

overlap.95 <- as.data.frame(intersect(distinct.spec$A, com.95.seqs$A))

# Comparison with P<0.5:
A <- substr(com.50$peptide, 12, nchar(com.50$peptide)-2) #extract characters from a cell with start and end point
B <- gsub("[[357.26]]", "", A)
B <- gsub("[[357.2]", "", B)
B <- gsub("1404]", "", B)
B <- gsub("1600", "", B)
B <- as.data.frame(B)
com.50.seqs <- as.data.frame(unique(B))
colnames(com.50.seqs) <- c('A')

overlap.50 <- as.data.frame(intersect(distinct.spec$A, com.50.seqs$A))







