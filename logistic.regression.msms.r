#------------------------------------------------------------------
#---------------------ADDING PEPs to msmsScan----------------------
#------------------------------------------------------------------

msms.scan <- read.delim('C:/Users/nslavov/Desktop/R-LogReg/msmsScans.txt', header=TRUE, stringsAsFactors=FALSE)
msms <- read.delim('C:/Users/nslavov/Desktop/R-LogReg/msms.txt', header=TRUE, stringsAsFactors=FALSE)

exps <- as.character(unique(msms.scan$Raw.file))

#msms.scan.pep <- data.frame(Raw.file=c(),Scan.number=c(),Retention.time=c(),Ion.injection.time=c(),Total.ion.current=c,Collision.energy=c(),Summations=c(),Base.peak.intensity=c(),Elapsed.time=c(),Identified=c(),MS.MS.IDs=c(),Sequence=c(),Length=c(),Filtered.peaks=c(),m.z=c(),Mass=c(),Charge=c(),Type=c(),Fragmentation=c(),Mass.analyzer=c(),Parent.intensity.fraction=c(),Fraction.of.total.spectrum=c(),Base.peak.fraction=c(),Precursor.full.scan.number=c(),Precursor.intensity=c(),Precursor.apex.fraction=c(),Precursor.apex.offset=c(),Precursor.apex.offset.time=c(),Scan.event.number=c(),Modifications=c(),Modified.sequence=c(),Proteins=c(),Score=c(),Intens.Comp.Factor=c(),CTCD.Comp=c(),RawOvFtT=c(),AGC.Fill=c(),Scan.index=c(),MS.scan.index=c(),MS.scan.number=c())

msms.scan.pep <- msms.scan[1,] #for having the headers
msms.scan.pep$PEP <- c(1) #for having PEP column
for (i in exps) {
  mssc.exp <- subset(msms.scan, Raw.file == i)
  ms.exp <- subset(msms, Raw.file == i)
  mssc.exp$PEP <- ms.exp$PEP[match(mssc.exp$Scan.index, ms.exp$Scan.index)]
  msms.scan.pep <- rbind(msms.scan.pep,mssc.exp)
}
msms.scan.pep <- msms.scan.pep[-1,] #deleting the first row

write.csv(msms.scan.pep, 'C:/Users/nslavov/Desktop/R-LogReg/msmsScan-With PEP.csv', row.names = FALSE)


#------------------------------------------------------------------
#------------------Creating PEP Categories-------------------------
#------------------------------------------------------------------

MSc <- read.csv('C:/Users/nslavov/Desktop/R-LogReg/msmsScan-With PEP.csv', header=TRUE, stringsAsFactors=FALSE)
msmsSc <- subset(MSc, Fragmentation=="HCD") #Keeping only HCDs
msmsSc <- subset(msmsSc,select=c('Retention.time', 'Base.peak.intensity', 'Base.peak.fraction', 'Ion.injection.time', 'Elapsed.time', 'm.z', 'Mass', 'Charge', 'Identified', 'PEP'))
msmsSc <- msmsSc[ ! msmsSc$Charge %in% c(0), ] #Remove scans with zero charges
msmsSc <- subset(msmsSc, PEP<=1) #Remove scans with PEP>1

for (i in 1:nrow(msmsSc)) {
  if (msmsSc$Identified[i]=="-"){msmsSc$PEP[i] <- 2}
}
write.csv(msmsSc, 'C:/Users/nslavov/Desktop/R-LogReg/msmsScan-With PEP and corrections.csv')

install.packages('Amelia')
require(Amelia)
missmap(msmsSc) # Check the map of missing values in the data frame


#------------CREATING CATEGORIES-----------------
for (i in 1:nrow(msmsSc)) {
  if (msmsSc$PEP[i]<=0.02){msmsSc$ID[i] <- 1}
  else if (msmsSc$PEP[i]>0.02 & msmsSc$PEP[i]<=0.05){msmsSc$ID[i] <- 2}
  else if (msmsSc$PEP[i]>0.05 & msmsSc$PEP[i]<=0.1){msmsSc$ID[i] <- 3}
  else if (msmsSc$PEP[i]>0.1 & msmsSc$PEP[i]<=1){msmsSc$ID[i] <- 4}
  else if (msmsSc$PEP[i]==2){msmsSc$ID[i] <- 0}
}
write.csv(msmsSc, 'C:/Users/nslavov/Desktop/R-LogReg/msmsScan-With PEP and corrections and ID column.csv')

#------------------------------------------------------------------
#----------------------Box Plots-----------------------------------
#------------------------------------------------------------------
msmsSc <- read.csv('C:/Users/nslavov/Desktop/R-LogReg/msmsScan-With PEP and corrections and ID column.csv', header=TRUE, stringsAsFactors=FALSE)

require(ggplot2)
#boxplot(Charge ~ iden, data = msmsSc, ylab = "Charge")
#boxplot(m.z ~ iden, data = msmsSc, ylab = "Charge")

msmsSc$ID <- as.factor(msmsSc$ID)

xlabs <- paste(levels(msmsSc$ID),"\n(N=",format(as.numeric(table(msmsSc$ID)), scientific = TRUE),")",sep="") #Scientific Notation

png('C:/Users/nslavov/Desktop/R-LogReg/charge-boxplot.png', width=8, height=6, units="in", res=300)
ggplot(msmsSc, aes(x=ID, y=Charge, fill=ID)) + 
  geom_boxplot()+xlab('Identification')+ylab('Charge')+
  theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+ #legend.position="none" removes the legend
  stat_summary(fun.y="mean", geom="point", size=2, position=position_dodge(width=0.75), color="white")+ #Adding the mean on the boxplots-White dot
  scale_x_discrete(labels=xlabs)+ #Adding the number of observations on the x axis
  scale_fill_discrete(labels=c("Unidentified", "[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
dev.off()

png('C:/Users/nslavov/Desktop/R-LogReg/precursor mz-boxplot.png', width=8, height=6, units="in", res=300)
ggplot(msmsSc, aes(x=ID, y=m.z, fill=ID)) + 
  geom_boxplot()+xlab('Identification')+ylab('Precursor m/z')+
  theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+ #legend.position="none" removes the legend
  stat_summary(fun.y="mean", geom="point", size=2, position=position_dodge(width=0.75), color="white")+ #Adding the mean on the boxplots-White dot
  scale_x_discrete(labels=xlabs)+ #Adding the number of observations on the x axis
  scale_fill_discrete(labels=c("Unidentified", "[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
dev.off()

png('C:/Users/nslavov/Desktop/R-LogReg/mass-boxplot.png', width=8, height=6, units="in", res=300)
ggplot(msmsSc, aes(x=ID, y=Mass, fill=ID)) + 
  geom_boxplot()+xlab('Identification')+ylab('Mass')+
  theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+ #legend.position="none" removes the legend
  stat_summary(fun.y="mean", geom="point", size=2, position=position_dodge(width=0.75), color="white")+ #Adding the mean on the boxplots-White dot
  scale_x_discrete(labels=xlabs)+ #Adding the number of observations on the x axis
  scale_fill_discrete(labels=c("Unidentified", "[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
dev.off()

png('C:/Users/nslavov/Desktop/R-LogReg/elapsed time-boxplot.png', width=8, height=6, units="in", res=300)
ggplot(msmsSc, aes(x=ID, y=Elapsed.time, fill=ID)) + 
  geom_boxplot()+xlab('Identification')+ylab('Elapsed Time')+
  theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+ #legend.position="none" removes the legend
  stat_summary(fun.y="mean", geom="point", size=2, position=position_dodge(width=0.75), color="white")+ #Adding the mean on the boxplots-White dot
  scale_x_discrete(labels=xlabs)+ #Adding the number of observations on the x axis
  scale_fill_discrete(labels=c("Unidentified", "[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
dev.off()

png('C:/Users/nslavov/Desktop/R-LogReg/retention time-boxplot.png', width=8, height=6, units="in", res=300)
ggplot(msmsSc, aes(x=ID, y=Retention.time, fill=ID)) + 
  geom_boxplot()+xlab('Identification')+ylab('Retention Time')+
  theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+ #legend.position="none" removes the legend
  stat_summary(fun.y="mean", geom="point", size=2, position=position_dodge(width=0.75), color="white")+ #Adding the mean on the boxplots-White dot
  scale_x_discrete(labels=xlabs)+ #Adding the number of observations on the x axis
  scale_fill_discrete(labels=c("Unidentified", "[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
dev.off()

png('C:/Users/nslavov/Desktop/R-LogReg/ion injection time-boxplot.png', width=8, height=6, units="in", res=300)
ggplot(msmsSc, aes(x=ID, y=Ion.injection.time, fill=ID)) + 
  geom_boxplot()+xlab('Identification')+ylab('Ion injection time')+
  theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+ #legend.position="none" removes the legend
  stat_summary(fun.y="mean", geom="point", size=2, position=position_dodge(width=0.75), color="white")+ #Adding the mean on the boxplots-White dot
  scale_x_discrete(labels=xlabs)+ #Adding the number of observations on the x axis
  scale_fill_discrete(labels=c("Unidentified", "[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
dev.off()


#------------------------------------------------------------------
#----------------------ALL PEPTIDES--------------------------------
#------------------------------------------------------------------

all.pep <- read.delim('C:/Users/nslavov/Desktop/R-LogReg/allPeptides.txt', header=TRUE, stringsAsFactors=FALSE)
all.pep <- all.pep[all.pep$Raw.file %in% c('071816L_19A_10ul_210min_100nl', '160629AEL160624PMRSAM00092_18C2', '160712A_NC_19D', '160720A_NC_19A_180min_100nl', '160728A_NC_set22A_180min_50ID_100NL', '160801A_NC_set19B_180min_50ID_100NL_from100ul', '160808L_NC_set#23A_180min', '160810A_NC_set25B_180min_100NL', '160818A_NC_set#19A_180min_200NL_60C+N=3', '160822A_NC_set#19B_180min_200NL_60C+N=3_500ms'), ] 
all.pep <- subset(all.pep,select=c('Raw.file', "MSMS.Scan.Numbers", 'Retention.length..FWHM.'))
all.pep <- all.pep[!(all.pep$MSMS.Scan.Numbers==""), ]
install.packages('splitstackshape')
require(splitstackshape)
all.pep <- cSplit(all.pep, "MSMS.Scan.Numbers", sep = ";", direction = "long")
write.csv(all.pep, 'C:/Users/nslavov/Desktop/R-LogReg/all.pep.filtered.csv', row.names = FALSE)

msms <- read.delim('C:/Users/nslavov/Desktop/R-LogReg/msms.txt', header=TRUE, stringsAsFactors=FALSE)
msms <- msms[msms$Raw.file %in% c('071816L_19A_10ul_210min_100nl', '160629AEL160624PMRSAM00092_18C2', '160712A_NC_19D', '160720A_NC_19A_180min_100nl', '160728A_NC_set22A_180min_50ID_100NL', '160801A_NC_set19B_180min_50ID_100NL_from100ul', '160808L_NC_set#23A_180min', '160810A_NC_set25B_180min_100NL', '160818A_NC_set#19A_180min_200NL_60C+N=3', '160822A_NC_set#19B_180min_200NL_60C+N=3_500ms'), ] 
msms <- subset(msms, Fragmentation=="HCD")
msms <- subset(msms,select=c('Raw.file', 'Scan.number', 'PEP'))
write.csv(msms, 'C:/Users/nslavov/Desktop/R-LogReg/msms.filtered.csv', row.names = FALSE)

#------------ADDING PEPs to ALL PEPTIDES------------

msms <- read.csv('C:/Users/nslavov/Desktop/R-LogReg/msms.filtered.csv', header=TRUE, stringsAsFactors=FALSE)
all.pep <- read.csv('C:/Users/nslavov/Desktop/R-LogReg/all.pep.filtered.csv', header=TRUE, stringsAsFactors=FALSE)
experiments <- as.character(c('071816L_19A_10ul_210min_100nl', '160629AEL160624PMRSAM00092_18C2', '160712A_NC_19D', '160720A_NC_19A_180min_100nl', '160728A_NC_set22A_180min_50ID_100NL', '160801A_NC_set19B_180min_50ID_100NL_from100ul', '160808L_NC_set#23A_180min', '160810A_NC_set25B_180min_100NL', '160818A_NC_set#19A_180min_200NL_60C+N=3', '160822A_NC_set#19B_180min_200NL_60C+N=3_500ms'))

all.pept.pep <- all.pep[1,]
all.pept.pep$PEP <- c(1) 
for (i in experiments) {
  msms.exp <- subset(msms, Raw.file == i)
  all.pep.exp <- subset(all.pep, Raw.file == i)
  all.pep.exp$PEP <- msms.exp$PEP[match(all.pep.exp$MSMS.Scan.Numbers, msms.exp$Scan.number)]
  all.pept.pep <- rbind(all.pept.pep,all.pep.exp)
}
all.pept.pep <- all.pept.pep[-1,]

missmap(all.pept.pep)
all.pept.pep <- all.pept.pep[!is.na(all.pept.pep$PEP), ]
write.csv(all.pept.pep, 'C:/Users/nslavov/Desktop/R-LogReg/all.pep.filtered+PEPs.csv', row.names = FALSE)

#------------CREATING CATEGORIES-----------------

all.pept.pep <- subset(all.pept.pep, PEP<=1)
for (i in 1:nrow(all.pept.pep)) {
  if (all.pept.pep$PEP[i]<=0.02){all.pept.pep$ID[i] <- 1}
  else if (all.pept.pep$PEP[i]>0.02 & all.pept.pep$PEP[i]<=0.05){all.pept.pep$ID[i] <- 2}
  else if (all.pept.pep$PEP[i]>0.05 & all.pept.pep$PEP[i]<=0.1){all.pept.pep$ID[i] <- 3}
  else if (all.pept.pep$PEP[i]>0.1 & all.pept.pep$PEP[i]<=1){all.pept.pep$ID[i] <- 4}
}
write.csv(all.pept.pep, 'C:/Users/nslavov/Desktop/R-LogReg/all.pept-With PEP and corrections and ID column.csv')

#----------------BOX PLOTS----------------------
all.pept.pep <- read.csv('C:/Users/nslavov/Desktop/R-LogReg/all.pept-With PEP and corrections and ID column.csv', header=TRUE, stringsAsFactors=FALSE)
experiments <- as.character(c('071816L_19A_10ul_210min_100nl', '160629AEL160624PMRSAM00092_18C2', '160712A_NC_19D', '160720A_NC_19A_180min_100nl', '160728A_NC_set22A_180min_50ID_100NL', '160801A_NC_set19B_180min_50ID_100NL_from100ul', '160808L_NC_set#23A_180min', '160810A_NC_set25B_180min_100NL', '160818A_NC_set#19A_180min_200NL_60C+N=3', '160822A_NC_set#19B_180min_200NL_60C+N=3_500ms'))
all.pept.pep$ID <- as.factor(all.pept.pep$ID)

for (i in experiments) {
  pep.exp <- subset(all.pept.pep, Raw.file == i)
  png.loc <- paste("C:/Users/nslavov/Desktop/R-LogReg/", i, ".png", sep = "")
  
  xlabs <- paste(levels(pep.exp$ID),"\n(N=",format(as.numeric(table(pep.exp$ID)), scientific = TRUE),")",sep="")
  
  png(png.loc, width=8, height=6, units="in", res=300)
  ggplot(pep.exp, aes(x=ID, y=Retention.length..FWHM., fill=ID)) + 
    geom_boxplot()+xlab('Identification')+ylab('FWHM')+ggtitle(i)+
    theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+ #legend.position="none" removes the legend
    stat_summary(fun.y="mean", geom="point", size=2, position=position_dodge(width=0.75), color="white")+ #Adding the mean on the boxplots-White dot
    scale_x_discrete(labels=xlabs)+ #Adding the number of observations on the x axis
    scale_fill_discrete(labels=c("[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
  dev.off()
}

  
  
  #--TEST-----------------------------------------------
  msmsSc.1 <- subset(msmsSc, PEP<=0.02 | PEP==2)
  for (i in 1:nrow(msmsSc.1)) {
    if (msmsSc.1$PEP[i]<1.5){msmsSc.1$ID[i] <- 1}
    else {msmsSc.1$ID[i] <- 0}
  }
  write.csv(msmsSc.1, 'C:/Users/nslavov/Desktop/R-LogReg/msmsSc.1.csv', row.names = FALSE)
  
  msmsSc.2 <- subset(msmsSc, PEP<=0.05 | PEP==2)
  for (i in 1:nrow(msmsSc.2)) {
    if (msmsSc.2$PEP[i]<1.5){msmsSc.2$ID[i] <- 2}
    else {msmsSc.2$ID[i] <- 0}
  }
  write.csv(msmsSc.2, 'C:/Users/nslavov/Desktop/R-LogReg/msmsSc.2.csv', row.names = FALSE)
  
  msmsSc.3 <- subset(msmsSc, PEP<=0.1 | PEP==2)
  for (i in 1:nrow(msmsSc.3)) {
    if (msmsSc.3$PEP[i]<1.5){msmsSc.3$ID[i] <- 3}
    else {msmsSc.3$ID[i] <- 0}
  }
  write.csv(msmsSc.3, 'C:/Users/nslavov/Desktop/R-LogReg/msmsSc.3.csv', row.names = FALSE)
  
  msmsSc.4 <- subset(msmsSc, PEP>=0.1)
  for (i in 1:nrow(msmsSc.4)) {
    if (msmsSc.4$PEP[i]<1.5){msmsSc.4$ID[i] <- 4}
    else {msmsSc.4$ID[i] <- 0}
  }
  write.csv(msmsSc.4, 'C:/Users/nslavov/Desktop/R-LogReg/msmsSc.4.csv', row.names = FALSE)
  #----
  #----
  B <- as.data.frame(A$A)
  B[is.na(B)] <- 100
  A$A <- B
  
  A <- data.frame(A=c(100,20000,3000,454546,544))
  for (i in 1:nrow(A)) {
    if (A$A[i]==2){print(i)}
  }
  
  
  
#------------------------------------------------------------------
#-------------Finding the correlations-----------------------------
#------------------------------------------------------------------
msmsSc.rand <- msmsSc[sample(nrow(msmsSc)),] #Randomizing the rows
msmsSc2$Charge <- as.numeric(msmsSc2$Charge)
msmsSc2$Identified <- NULL
msmsSc2 <- msmsSc2[sample(nrow(msmsSc2)),]

model <- glm(iden ~ Base.peak.intensity + Elapsed.time + m.z + Mass, data = msmsSc2, family=binomial(link='logit'))
table(msmsSc2$iden,predict(model,type='response')>=0.5) #contingency table

Retention.time	<- cor(msmsSc$Retention.time, msmsSc$iden, method = "pearson")











