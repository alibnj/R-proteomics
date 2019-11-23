#-------------------------
#------SpectraST/Comet-------
#-------------------------

## Reading and cleaning
spec <- read.csv('C:/Inetpub/wwwroot/ISB/data/human.spec-v3.quantitation/new/160831A_NC_set#27X_1ul_180min_200NL_60C+N=3_250ms.csv', header=TRUE, stringsAsFactors=FALSE)
spec <- subset(spec, spec$probability>=0.99)
spec <- spec[c('peptide', 'libra1', 'libra2', 'libra3', 'libra4', 'libra5', 'libra6', 'libra7', 'libra9')]
spec <- spec[!(spec$libra1==0 & spec$libra2==0 & spec$libra3==0 & spec$libra4==0 & spec$libra5==0 & spec$libra6==0 & spec$libra7==0 & spec$libra9==0), ]

## Normalizing quantitations
spec$libra1 <- spec$libra1/mean(spec$libra1)
spec$libra2 <- spec$libra2/mean(spec$libra2)
spec$libra3 <- spec$libra3/mean(spec$libra3)
spec$libra4 <- spec$libra4/mean(spec$libra4)
spec$libra5 <- spec$libra5/mean(spec$libra5)
spec$libra6 <- spec$libra6/mean(spec$libra6)
spec$libra7 <- spec$libra7/mean(spec$libra7)
spec$libra9 <- spec$libra9/mean(spec$libra9)

#for (i in 1:nrow(spec)) {
#  spec[i,2:9] <- spec[i,2:9]/(mean(as.numeric(spec[i,2:9]))) #Normalizing the quants by deviding channel's intensities of each row to the mean of the intensities of that row
#}

## Collapsing the quantitations for each peptide
spec.uniq.quant <- data.frame(peptide = character(), ch1 = numeric(), ch2 = numeric(), ch3 = numeric(), ch4 = numeric(), ch5 = numeric(), ch6 = numeric(), ch7 = numeric(), ch9 = numeric(), stringsAsFactors=FALSE)
A <- data.frame(peptide = character(), ch1 = numeric(), ch2 = numeric(), ch3 = numeric(), ch4 = numeric(), ch5 = numeric(), ch6 = numeric(), ch7 = numeric(), ch9 = numeric(), stringsAsFactors=FALSE)
spec.peps <- as.character(unique(spec$peptide))
for (i in spec.peps) {
  spec.p <- subset(spec, spec$peptide==i)
  A[1,1] <- i
  A$ch1 <- median(spec.p$libra1)
  A$ch2 <- median(spec.p$libra2)
  A$ch3 <- median(spec.p$libra3)
  A$ch4 <- median(spec.p$libra4)
  A$ch5 <- median(spec.p$libra5)
  A$ch6 <- median(spec.p$libra6)
  A$ch7 <- median(spec.p$libra7)
  A$ch9 <- median(spec.p$libra9)
  spec.uniq.quant <- rbind(spec.uniq.quant, A)
}
#write.csv(spec.uniq.quant, 'C:/Inetpub/wwwroot/ISB/data/human.spec-v3.quantitation/quant.correlation/corr/all-0.99.uniq.quant-spec.csv', row.names = FALSE)




#-------------------------
#------Max Quant------------
#-------------------------

## Reading and Cleaning
mq <- read.delim('D:/Ali/SingleCells/RAW/combined/txt/evidence.txt', header=TRUE, stringsAsFactors=FALSE)
mq <- subset(mq, mq$PEP<=0.01)
#mq <- mq[c("Raw.file", "Sequence","Reporter.intensity.corrected.0","Reporter.intensity.corrected.1","Reporter.intensity.corrected.2","Reporter.intensity.corrected.3","Reporter.intensity.corrected.4","Reporter.intensity.corrected.5","Reporter.intensity.corrected.6","Reporter.intensity.corrected.8")]
mq <- mq[c("Raw.file", "Sequence","Reporter.intensity.0","Reporter.intensity.1","Reporter.intensity.2","Reporter.intensity.3","Reporter.intensity.4","Reporter.intensity.5","Reporter.intensity.6","Reporter.intensity.8")]
mq <- mq[!(mq$Reporter.intensity.0==0 & mq$Reporter.intensity.1==0 & mq$Reporter.intensity.2==0 & mq$Reporter.intensity.3==0 & mq$Reporter.intensity.4==0 & mq$Reporter.intensity.5==0 & mq$Reporter.intensity.6==0 & mq$Reporter.intensity.8==0), ]

mq <- subset(mq, mq$Raw.file=='160831A_NC_set#27X_1ul_180min_200NL_60C+N=3_250ms')
mq.peps <- unique(mq$Sequence)

## Normalizing Quantitations
#mq$Reporter.intensity.corrected.0 <- mq$Reporter.intensity.corrected.0/mean(mq$Reporter.intensity.corrected.0)
#mq$Reporter.intensity.corrected.1 <- mq$Reporter.intensity.corrected.1/mean(mq$Reporter.intensity.corrected.1)
#mq$Reporter.intensity.corrected.2 <- mq$Reporter.intensity.corrected.2/mean(mq$Reporter.intensity.corrected.2)
#mq$Reporter.intensity.corrected.3 <- mq$Reporter.intensity.corrected.3/mean(mq$Reporter.intensity.corrected.3)
#mq$Reporter.intensity.corrected.4 <- mq$Reporter.intensity.corrected.4/mean(mq$Reporter.intensity.corrected.4)
#mq$Reporter.intensity.corrected.5 <- mq$Reporter.intensity.corrected.5/mean(mq$Reporter.intensity.corrected.5)
#mq$Reporter.intensity.corrected.6 <- mq$Reporter.intensity.corrected.6/mean(mq$Reporter.intensity.corrected.6)
#mq$Reporter.intensity.corrected.8 <- mq$Reporter.intensity.corrected.8/mean(mq$Reporter.intensity.corrected.8)

mq$Reporter.intensity.0 <- mq$Reporter.intensity.0/mean(mq$Reporter.intensity.0)
mq$Reporter.intensity.1 <- mq$Reporter.intensity.1/mean(mq$Reporter.intensity.1)
mq$Reporter.intensity.2 <- mq$Reporter.intensity.2/mean(mq$Reporter.intensity.2)
mq$Reporter.intensity.3 <- mq$Reporter.intensity.3/mean(mq$Reporter.intensity.3)
mq$Reporter.intensity.4 <- mq$Reporter.intensity.4/mean(mq$Reporter.intensity.4)
mq$Reporter.intensity.5 <- mq$Reporter.intensity.5/mean(mq$Reporter.intensity.5)
mq$Reporter.intensity.6 <- mq$Reporter.intensity.6/mean(mq$Reporter.intensity.6)
mq$Reporter.intensity.8 <- mq$Reporter.intensity.8/mean(mq$Reporter.intensity.8)

## Collapsing the quantitations for each peptide
mq.uniq.quant <- data.frame(peptide = character(), ch1 = numeric(), ch2 = numeric(), ch3 = numeric(), ch4 = numeric(), ch5 = numeric(), ch6 = numeric(), ch7 = numeric(), ch9 = numeric(), stringsAsFactors=FALSE)
A <- data.frame(peptide = character(), ch1 = numeric(), ch2 = numeric(), ch3 = numeric(), ch4 = numeric(), ch5 = numeric(), ch6 = numeric(), ch7 = numeric(), ch9 = numeric(), stringsAsFactors=FALSE)
mq.peps <- as.character(unique(mq$Sequence))
#for (i in mq.peps) {
#  mq.p <- subset(mq, mq$Sequence==i)
#  A[1, 1] <- i
#  A$ch1 <- median(mq.p$Reporter.intensity.corrected.0)
#  A$ch2 <- median(mq.p$Reporter.intensity.corrected.1)
#  A$ch3 <- median(mq.p$Reporter.intensity.corrected.2)
#  A$ch4 <- median(mq.p$Reporter.intensity.corrected.3)
#  A$ch5 <- median(mq.p$Reporter.intensity.corrected.4)
#  A$ch6 <- median(mq.p$Reporter.intensity.corrected.5)
#  A$ch7 <- median(mq.p$Reporter.intensity.corrected.6)
#  A$ch9 <- median(mq.p$Reporter.intensity.corrected.8)
#  mq.uniq.quant <- rbind(mq.uniq.quant, A)
#}

for (i in mq.peps) {
  mq.p <- subset(mq, mq$Sequence==i)
  A[1, 1] <- i
  A$ch1 <- median(mq.p$Reporter.intensity.0)
  A$ch2 <- median(mq.p$Reporter.intensity.1)
  A$ch3 <- median(mq.p$Reporter.intensity.2)
  A$ch4 <- median(mq.p$Reporter.intensity.3)
  A$ch5 <- median(mq.p$Reporter.intensity.4)
  A$ch6 <- median(mq.p$Reporter.intensity.5)
  A$ch7 <- median(mq.p$Reporter.intensity.6)
  A$ch9 <- median(mq.p$Reporter.intensity.8)
  mq.uniq.quant <- rbind(mq.uniq.quant, A)
}

mq.uniq.quant.spec <- mq.uniq.quant[mq.uniq.quant$peptide %in% spec.peps, ] #Keeping the overlap with spec
mq.sp.peps <- unique(mq.uniq.quant.spec$peptide)
spec.uniq.quant.mq <- spec.uniq.quant[spec.uniq.quant$peptide %in% mq.sp.peps, ] #Keeping the shared peptides with mq

#write.csv(mq.uniq.quant, 'C:/Inetpub/wwwroot/ISB/data/human.spec-v3.quantitation/quant.correlation/corr/all-correct.uniq.quant-mq-only.spec.overlaps.csv', row.names = FALSE)
#write.csv(spec.uniq.quant.mq, 'C:/Inetpub/wwwroot/ISB/data/human.spec-v3.quantitation/quant.correlation/corr/all-correct.uniq.quant-spec-only.mq.overlaps.csv', row.names = FALSE)


#-------------------------
#------Comet----------------
#-------------------------
## Reading and cleaning
com <- read.csv('C:/Inetpub/wwwroot/ISB/data/human.comet-v2.quantitation/new/160831A_NC_set#19A_180min_200NL_60C+N=3_250ms.csv', header=TRUE, stringsAsFactors=FALSE)
com <- subset(com, com$probability>=0.99)
com <- com[c('peptide', 'libra1', 'libra2', 'libra3', 'libra4', 'libra5', 'libra6', 'libra7', 'libra9')]
com <- com[!(com$libra1==0 & com$libra2==0 & com$libra3==0 & com$libra4==0 & com$libra5==0 & com$libra6==0 & com$libra7==0 & com$libra9==0), ]

## Normalizing quantitations
com$libra1 <- com$libra1/mean(com$libra1)
com$libra2 <- com$libra2/mean(com$libra2)
com$libra3 <- com$libra3/mean(com$libra3)
com$libra4 <- com$libra4/mean(com$libra4)
com$libra5 <- com$libra5/mean(com$libra5)
com$libra6 <- com$libra6/mean(com$libra6)
com$libra7 <- com$libra7/mean(com$libra7)
com$libra9 <- com$libra9/mean(com$libra9)

## Collapsing the quantitations for each peptide
com.uniq.quant <- data.frame(peptide = character(), ch1 = numeric(), ch2 = numeric(), ch3 = numeric(), ch4 = numeric(), ch5 = numeric(), ch6 = numeric(), ch7 = numeric(), ch9 = numeric(), stringsAsFactors=FALSE)
A <- data.frame(peptide = character(), ch1 = numeric(), ch2 = numeric(), ch3 = numeric(), ch4 = numeric(), ch5 = numeric(), ch6 = numeric(), ch7 = numeric(), ch9 = numeric(), stringsAsFactors=FALSE)
com.peps <- as.character(unique(com$peptide))
for (i in com.peps) {
  com.p <- subset(com, com$peptide==i)
  A[1,1] <- i
  A$ch1 <- median(com.p$libra1)
  A$ch2 <- median(com.p$libra2)
  A$ch3 <- median(com.p$libra3)
  A$ch4 <- median(com.p$libra4)
  A$ch5 <- median(com.p$libra5)
  A$ch6 <- median(com.p$libra6)
  A$ch7 <- median(com.p$libra7)
  A$ch9 <- median(com.p$libra9)
  com.uniq.quant <- rbind(com.uniq.quant, A)
}
write.csv(com.uniq.quant, 'C:/Inetpub/wwwroot/ISB/data/human.com-v3.quantitation/quant.correlation/corr/all-0.99.uniq.quant-com.csv', row.names = FALSE)

com.uniq.quant.spec <- com.uniq.quant[com.uniq.quant$peptide %in% spec.peps, ] #Keeping the overlap with spec
com.sp.peps <- unique(com.uniq.quant.spec$peptide)
spec.uniq.quant.com <- spec.uniq.quant[spec.uniq.quant$peptide %in% com.sp.peps, ]


#-------------------------
#------Correlations---------
#-------------------------

## MQ(PEP=0.02) vs. SpectraST (Libra with contributions+Mass tol=0.01+Min P=0.8)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Sorting both data frames to compare same peptides row-wise:
mq.q <- mq.uniq.quant.spec[order(mq.uniq.quant.spec$peptide), ]
rownames(mq.q) <- NULL #resetting rows index
mq.q$peptide <- NULL

mq.q <- as.data.frame(mq.q)
mq.q$ch1 <- mq.q$ch1/sum(mq.q$ch1)
mq.q$ch2 <- mq.q$ch2/sum(mq.q$ch2)
mq.q$ch3 <- mq.q$ch3/sum(mq.q$ch3)
mq.q$ch4 <- mq.q$ch4/sum(mq.q$ch4)
mq.q$ch5 <- mq.q$ch5/sum(mq.q$ch5)
mq.q$ch6 <- mq.q$ch6/sum(mq.q$ch6)
mq.q$ch7 <- mq.q$ch7/sum(mq.q$ch7)
mq.q$ch9 <- mq.q$ch9/sum(mq.q$ch9)

mq.q$ch1 <- scale(mq.q$ch1)
mq.q$ch2 <- scale(mq.q$ch2)
mq.q$ch3 <- scale(mq.q$ch3)
mq.q$ch4 <- scale(mq.q$ch4)
mq.q$ch5 <- scale(mq.q$ch5)
mq.q$ch6 <- scale(mq.q$ch6)
mq.q$ch7 <- scale(mq.q$ch7)
mq.q$ch9 <- scale(mq.q$ch9)

mq.q$ch1 <- mq.q$ch1/mean(mq.q$ch9)
mq.q$ch2 <- mq.q$ch2/mean(mq.q$ch9)
mq.q$ch3 <- mq.q$ch3/mean(mq.q$ch9)
mq.q$ch4 <- mq.q$ch4/mean(mq.q$ch9)
mq.q$ch5 <- mq.q$ch5/mean(mq.q$ch9)
mq.q$ch6 <- mq.q$ch6/mean(mq.q$ch9)
mq.q$ch7 <- mq.q$ch7/mean(mq.q$ch9)
mq.q$ch9 <- mq.q$ch9/mean(mq.q$ch9)

mq.q$ch1 <- (mq.q$ch1-min(mq.q$ch1)/max(mq.q$ch1)-min(mq.q$ch1))
mq.q$ch2 <- (mq.q$ch2-min(mq.q$ch2)/max(mq.q$ch2)-min(mq.q$ch2))
mq.q$ch3 <- (mq.q$ch3-min(mq.q$ch3)/max(mq.q$ch3)-min(mq.q$ch3))
mq.q$ch4 <- (mq.q$ch4-min(mq.q$ch4)/max(mq.q$ch4)-min(mq.q$ch4))
mq.q$ch5 <- (mq.q$ch5-min(mq.q$ch5)/max(mq.q$ch5)-min(mq.q$ch5))
mq.q$ch6 <- (mq.q$ch6-min(mq.q$ch6)/max(mq.q$ch6)-min(mq.q$ch6))
mq.q$ch7 <- (mq.q$ch7-min(mq.q$ch7)/max(mq.q$ch7)-min(mq.q$ch7))
mq.q$ch9 <- (mq.q$ch9-min(mq.q$ch9)/max(mq.q$ch9)-min(mq.q$ch9))


for (i in 1:nrow(mq.q)) {
  mq.q[i,] <- mq.q[i,]/(sum(as.numeric(mq.q[i,]))) #Normalizing the quants by deviding channel's intensities of each row to the mean of the intensities of that row
}
#----------------------------------
###--------LIBRA------------------->>>>
#----------------------------------
m.mq.ch1 <- mean(mq.q$ch1)
m.mq.ch2 <- mean(mq.q$ch2)
m.mq.ch3 <- mean(mq.q$ch3)
m.mq.ch4 <- mean(mq.q$ch4)
m.mq.ch5 <- mean(mq.q$ch5)
m.mq.ch6 <- mean(mq.q$ch6)
m.mq.ch7 <- mean(mq.q$ch7)
m.mq.ch9 <- mean(mq.q$ch9)

sd.mq.ch1 <- sd(mq.q$ch1)
sd.mq.ch2 <- sd(mq.q$ch2)
sd.mq.ch3 <- sd(mq.q$ch3)
sd.mq.ch4 <- sd(mq.q$ch4)
sd.mq.ch5 <- sd(mq.q$ch5)
sd.mq.ch6 <- sd(mq.q$ch6)
sd.mq.ch7 <- sd(mq.q$ch7)
sd.mq.ch9 <- sd(mq.q$ch9)

for (i in 1:nrow(mq.q)) {
  if (mq.q$ch1[i]>((m.mq.ch1)+2*sd(mq.q$ch1)) | mq.q$ch1[i]<((m.mq.ch1)-2*sd(mq.q$ch1))){mq.q$ch1[i] <- 0}
  if (mq.q$ch2[i]>((m.mq.ch2)+2*sd(mq.q$ch2)) | mq.q$ch2[i]<((m.mq.ch2)-2*sd(mq.q$ch2))){mq.q$ch2[i] <- 0}
  if (mq.q$ch3[i]>((m.mq.ch3)+2*sd(mq.q$ch3)) | mq.q$ch3[i]<((m.mq.ch3)-2*sd(mq.q$ch3))){mq.q$ch3[i] <- 0}
  if (mq.q$ch4[i]>((m.mq.ch4)+2*sd(mq.q$ch4)) | mq.q$ch4[i]<((m.mq.ch4)-2*sd(mq.q$ch4))){mq.q$ch4[i] <- 0}
  if (mq.q$ch5[i]>((m.mq.ch5)+2*sd(mq.q$ch5)) | mq.q$ch5[i]<((m.mq.ch5)-2*sd(mq.q$ch5))){mq.q$ch5[i] <- 0}
  if (mq.q$ch6[i]>((m.mq.ch6)+2*sd(mq.q$ch6)) | mq.q$ch6[i]<((m.mq.ch6)-2*sd(mq.q$ch6))){mq.q$ch6[i] <- 0}
  if (mq.q$ch7[i]>((m.mq.ch7)+2*sd(mq.q$ch7)) | mq.q$ch7[i]<((m.mq.ch7)-2*sd(mq.q$ch7))){mq.q$ch7[i] <- 0}
  if (mq.q$ch9[i]>((m.mq.ch9)+2*sd(mq.q$ch9)) | mq.q$ch9[i]<((m.mq.ch9)-2*sd(mq.q$ch9))){mq.q$ch9[i] <- 0}
}

#-----------------------------------------------------------------------------------------------------

spec.q <- spec.uniq.quant.mq[order(spec.uniq.quant.mq$peptide), ]
rownames(spec.q) <- NULL
spec.q$peptide <- NULL

spec.q <- as.data.frame(spec.q)
spec.q$ch1 <- spec.q$ch1/sum(spec.q$ch1)
spec.q$ch2 <- spec.q$ch2/sum(spec.q$ch2)
spec.q$ch3 <- spec.q$ch3/sum(spec.q$ch3)
spec.q$ch4 <- spec.q$ch4/sum(spec.q$ch4)
spec.q$ch5 <- spec.q$ch5/sum(spec.q$ch5)
spec.q$ch6 <- spec.q$ch6/sum(spec.q$ch6)
spec.q$ch7 <- spec.q$ch7/sum(spec.q$ch7)
spec.q$ch9 <- spec.q$ch9/sum(spec.q$ch9)

spec.q$ch1 <- scale(spec.q$ch1)
spec.q$ch2 <- scale(spec.q$ch2)
spec.q$ch3 <- scale(spec.q$ch3)
spec.q$ch4 <- scale(spec.q$ch4)
spec.q$ch5 <- scale(spec.q$ch5)
spec.q$ch6 <- scale(spec.q$ch6)
spec.q$ch7 <- scale(spec.q$ch7)
spec.q$ch9 <- scale(spec.q$ch9)

spec.q$ch1 <- spec.q$ch1/mean(spec.q$ch9)
spec.q$ch2 <- spec.q$ch2/mean(spec.q$ch9)
spec.q$ch3 <- spec.q$ch3/mean(spec.q$ch9)
spec.q$ch4 <- spec.q$ch4/mean(spec.q$ch9)
spec.q$ch5 <- spec.q$ch5/mean(spec.q$ch9)
spec.q$ch6 <- spec.q$ch6/mean(spec.q$ch9)
spec.q$ch7 <- spec.q$ch7/mean(spec.q$ch9)
spec.q$ch9 <- spec.q$ch9/mean(spec.q$ch9)

spec.q$ch1 <- (spec.q$ch1-min(spec.q$ch1)/max(spec.q$ch1)-min(spec.q$ch1))
spec.q$ch2 <- (spec.q$ch2-min(spec.q$ch2)/max(spec.q$ch2)-min(spec.q$ch2))
spec.q$ch3 <- (spec.q$ch3-min(spec.q$ch3)/max(spec.q$ch3)-min(spec.q$ch3))
spec.q$ch4 <- (spec.q$ch4-min(spec.q$ch4)/max(spec.q$ch4)-min(spec.q$ch4))
spec.q$ch5 <- (spec.q$ch5-min(spec.q$ch5)/max(spec.q$ch5)-min(spec.q$ch5))
spec.q$ch6 <- (spec.q$ch6-min(spec.q$ch6)/max(spec.q$ch6)-min(spec.q$ch6))
spec.q$ch7 <- (spec.q$ch7-min(spec.q$ch7)/max(spec.q$ch7)-min(spec.q$ch7))
spec.q$ch9 <- (spec.q$ch9-min(spec.q$ch9)/max(spec.q$ch9)-min(spec.q$ch9))

for (i in 1:nrow(spec.q)) {
  spec.q[i,] <- spec.q[i,]/(sum(as.numeric(spec.q[i,]))) #Normalizing the quants by deviding channel's intensities of each row to the mean of the intensities of that row
}

spec.q <- scale(spec.q)
mq.q <- scale(mq.q)

spec.q <- (spec.q-min(spec.q))/(max(spec.q)-min(spec.q))
mq.q <- (mq.q-min(mq.q))/(max(mq.q)-min(mq.q))

#----------------------------------
###--------LIBRA------------------->>>>
#----------------------------------
m.sp.ch1 <- mean(spec.q$ch1)
m.sp.ch2 <- mean(spec.q$ch2)
m.sp.ch3 <- mean(spec.q$ch3)
m.sp.ch4 <- mean(spec.q$ch4)
m.sp.ch5 <- mean(spec.q$ch5)
m.sp.ch6 <- mean(spec.q$ch6)
m.sp.ch7 <- mean(spec.q$ch7)
m.sp.ch9 <- mean(spec.q$ch9)

sd.sp.ch1 <- sd(spec.q$ch1)
sd.sp.ch2 <- sd(spec.q$ch2)
sd.sp.ch3 <- sd(spec.q$ch3)
sd.sp.ch4 <- sd(spec.q$ch4)
sd.sp.ch5 <- sd(spec.q$ch5)
sd.sp.ch6 <- sd(spec.q$ch6)
sd.sp.ch7 <- sd(spec.q$ch7)
sd.sp.ch9 <- sd(spec.q$ch9)

for (i in 1:nrow(spec.q)) {
  if (spec.q$ch1[i]>((m.sp.ch1)+2*sd(spec.q$ch1)) | spec.q$ch1[i]<((m.sp.ch1)-2*sd(spec.q$ch1))){spec.q$ch1[i] <- 0}
  if (spec.q$ch2[i]>((m.sp.ch2)+2*sd(spec.q$ch2)) | spec.q$ch2[i]<((m.sp.ch2)-2*sd(spec.q$ch2))){spec.q$ch2[i] <- 0}
  if (spec.q$ch3[i]>((m.sp.ch3)+2*sd(spec.q$ch3)) | spec.q$ch3[i]<((m.sp.ch3)-2*sd(spec.q$ch3))){spec.q$ch3[i] <- 0}
  if (spec.q$ch4[i]>((m.sp.ch4)+2*sd(spec.q$ch4)) | spec.q$ch4[i]<((m.sp.ch4)-2*sd(spec.q$ch4))){spec.q$ch4[i] <- 0}
  if (spec.q$ch5[i]>((m.sp.ch5)+2*sd(spec.q$ch5)) | spec.q$ch5[i]<((m.sp.ch5)-2*sd(spec.q$ch5))){spec.q$ch5[i] <- 0}
  if (spec.q$ch6[i]>((m.sp.ch6)+2*sd(spec.q$ch6)) | spec.q$ch6[i]<((m.sp.ch6)-2*sd(spec.q$ch6))){spec.q$ch6[i] <- 0}
  if (spec.q$ch7[i]>((m.sp.ch7)+2*sd(spec.q$ch7)) | spec.q$ch7[i]<((m.sp.ch7)-2*sd(spec.q$ch7))){spec.q$ch7[i] <- 0}
  if (spec.q$ch9[i]>((m.sp.ch9)+2*sd(spec.q$ch9)) | spec.q$ch9[i]<((m.sp.ch9)-2*sd(spec.q$ch9))){spec.q$ch9[i] <- 0}
}


#----------------------------------
# Finding the corelation:
#----------------------------------
## SpectraST vs. Max Quant
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
mq.q <- as.matrix(mq.q)
spec.q <- as.matrix(spec.q)
Z <- sapply(seq.int(dim(mq.q)[1]), function(i) cor(mq.q[i,], spec.q[i,], method = c("pearson"))) #Comparing the correlations row by row
hist(Z, breaks=seq(-1,1,by=0.1), xlab = 'Correlation between MQ & SpectraST', ylab = '# Unique peptides', main = NULL)
mean(Z)
Z <- cor(as.numeric(mq.q[,2]), as.numeric(spec.q[,2]))
mq.q <- mq.q %*% t(solve(corr.m))
spec.q <- spec.q %*% t(solve(corr.m))


## Comet vs. SpectraST
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Sorting both data frames to compare same peptides row-wise:
com.q <- com.uniq.quant.spec[order(com.uniq.quant.spec$peptide), ]
rownames(com.q) <- NULL #resetting rows index
com.q$peptide <- NULL

spec.q <- spec.uniq.quant.com[order(spec.uniq.quant.com$peptide), ]
rownames(spec.q) <- NULL
spec.q$peptide <- NULL

# Finding the corelation:
com.q <- as.matrix(com.q)
spec.q <- as.matrix(spec.q)
Z <- sapply(seq.int(dim(com.q)[1]), function(i) cor(com.q[i,], spec.q[i,])) #Comparing the correlations row by row
hist(Z, breaks=seq(-1,1,by=0.1), xlab = 'Correlation between Comet & SpectraST', ylab = '# Unique peptides', main = NULL)
mean(Z)


## Comet vs. MQ
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
com.uniq.quant.mq <- com.uniq.quant[com.uniq.quant$peptide %in% mq.peps, ] #Keeping the overlap with spec
com.sp.peps <- unique(com.uniq.quant.mq$peptide)
mq.uniq.quant.com <- mq.uniq.quant[mq.uniq.quant$peptide %in% com.sp.peps, ]


# Sorting both data frames to compare same peptides row-wise:
com.q <- com.uniq.quant.mq[order(com.uniq.quant.mq$peptide), ]
rownames(com.q) <- NULL #resetting rows index
com.q$peptide <- NULL

mq.q <- mq.uniq.quant.com[order(mq.uniq.quant.com$peptide), ]
rownames(mq.q) <- NULL #resetting rows index
mq.q$peptide <- NULL

# Finding the corelation:
com.q <- as.matrix(com.q)
mq.q <- as.matrix(mq.q)
Z <- sapply(seq.int(dim(com.q)[1]), function(i) cor(com.q[i,], mq.q[i,])) #Comparing the correlations row by row
hist(Z, breaks=seq(-1,1,by=0.1), xlab = 'Correlation between Comet & MaxQuant', ylab = '# Unique peptides', main = NULL)
mean(Z)

