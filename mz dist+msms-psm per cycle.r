MSc <- read.csv('C:/Users/nslavov/Desktop/R-LogReg/msmsScan-With PEP.csv', header=TRUE, stringsAsFactors=FALSE)
msmsSc <- subset(MSc, Fragmentation=="HCD")
msmsSc <- subset(msmsSc,select=c('MS.scan.index', 'Retention.time', 'Base.peak.intensity', 'Base.peak.fraction', 'Ion.injection.time', 'Elapsed.time', 'm.z', 'Mass', 'Charge', 'Identified', 'PEP'))
#msmsSc <- msmsSc[ ! msmsSc$Charge %in% c(0), ] #Remove scans with zero charges
#msmsSc <- subset(msmsSc, PEP<=1) #Remove scans with PEP>1

msmsSc[msmsSc$Identified=="-","PEP"]<-2
msmsSc[msmsSc$Charge==0,"PEP"]<-2
msmsSc[is.na(msmsSc$PEP),"PEP"]<- 500
msmsSc[msmsSc$PEP>1,"PEP"]<-2

A <- data.frame(A=c(1,2,3,4,5,6,7,8,9,10))
A[A$A>1 & A$A<5,"B"]<-'Z'

#This loops takes for everto finish!!!
#for (i in 1:nrow(msmsSc)) {
#  if (msmsSc$Identified[i]=="-" | msmsSc$Charge==0 | msmsSc$PEP>1){msmsSc$PEP[i] <- 2}
#}

write.csv(msmsSc, 'C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/msmsScan-With PEP and (-&charge0&PEP above 1) uniden.csv', row.names = FALSE)

install.packages('Amelia')
require(Amelia)
missmap(msmsSc)

msmsSc[msmsSc$PEP==2,"ID"] <- 0
msmsSc[msmsSc$PEP<=0.02,"ID"] <- 1
msmsSc[msmsSc$PEP>0.02 & msmsSc$PEP<=0.05,"ID"] <- 2
msmsSc[msmsSc$PEP>0.05 & msmsSc$PEP<=0.1,"ID"] <- 3
msmsSc[msmsSc$PEP>0.1 & msmsSc$PEP<=1,"ID"] <- 4

write.csv(msmsSc, 'C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/msmsScan-With PEP and corrections and ID column.csv', row.names = FALSE)


msmsSc <- read.csv('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/msmsScan-With PEP and corrections and ID column.csv', header=TRUE, stringsAsFactors=FALSE)


cat0 <- subset(msmsSc, ID=='0')
cat0.cycle <- unique(cat0$MS.scan.index)
A <- data.frame(MS.scan.index=numeric(), Freq=numeric())
C0 <- data.frame(MS.scan.index=numeric(), Freq=numeric())
for (i in cat0.cycle){
  AA <- subset(cat0, MS.scan.index==i)
  A[1,1] <- i
  A[1,2] <- nrow(AA)
  C0 <- rbind(C0, A)
}
C0$ID <- 'Unidentified'
write.csv(C0, 'C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C0.csv', row.names = FALSE)



cat1 <- subset(msmsSc, ID=='1')
cat1.cycle <- unique(cat1$MS.scan.index)
A <- data.frame(MS.scan.index=numeric(), Freq=numeric())
C1 <- data.frame(MS.scan.index=numeric(), Freq=numeric())
for (i in cat1.cycle){
  AA <- subset(cat1, MS.scan.index==i)
  A[1,1] <- i
  A[1,2] <- nrow(AA)
  C1 <- rbind(C1, A)
}
C1$ID <- '[0, 0.02]'
write.csv(C1, 'C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C1.csv', row.names = FALSE)

cat2 <- subset(msmsSc, ID=='2')
cat2.cycle <- unique(cat2$MS.scan.index)
A <- data.frame(MS.scan.index=numeric(), Freq=numeric())
C2 <- data.frame(MS.scan.index=numeric(), Freq=numeric())
for (i in cat2.cycle){
  AA <- subset(cat2, MS.scan.index==i)
  A[1,1] <- i
  A[1,2] <- nrow(AA)
  C2 <- rbind(C2, A)
}
C2$ID <- '(0.02, 0.05]'
write.csv(C2, 'C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C2.csv', row.names = FALSE)

cat3 <- subset(msmsSc, ID=='3')
cat3.cycle <- unique(cat3$MS.scan.index)
A <- data.frame(MS.scan.index=numeric(), Freq=numeric())
C3 <- data.frame(MS.scan.index=numeric(), Freq=numeric())
for (i in cat3.cycle){
  AA <- subset(cat3, MS.scan.index==i)
  A[1,1] <- i
  A[1,2] <- nrow(AA)
  C3 <- rbind(C3, A)
}
C3$ID <- '(0.05, 0.1]'
write.csv(C3, 'C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C3.csv', row.names = FALSE)

cat4 <- subset(msmsSc, ID=='4')
cat4.cycle <- unique(cat4$MS.scan.index)
A <- data.frame(MS.scan.index=numeric(), Freq=numeric())
C4 <- data.frame(MS.scan.index=numeric(), Freq=numeric())
for (i in cat4.cycle){
  AA <- subset(cat4, MS.scan.index==i)
  A[1,1] <- i
  A[1,2] <- nrow(AA)
  C4 <- rbind(C4, A)
}
C4$ID <- '(0.1, 1]'
write.csv(C4, 'C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C4.csv', row.names = FALSE)


msmsSc.cycle <- unique(msmsSc$MS.scan.index)
A <- data.frame(MS.scan.index=numeric(), Freq=numeric())
C5 <- data.frame(MS.scan.index=numeric(), Freq=numeric())
for (i in msmsSc.cycle){
  AA <- subset(msmsSc, MS.scan.index==i)
  A[1,1] <- i
  A[1,2] <- nrow(AA)
  C5 <- rbind(C5, A)
}
C5$ID <- 'All MSMS'
write.csv(C5, 'C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C5.csv', row.names = FALSE)

C0 <- read.csv('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C0.csv', header=TRUE, stringsAsFactors=FALSE)
C1 <- read.csv('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C1.csv', header=TRUE, stringsAsFactors=FALSE)
C2 <- read.csv('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C2.csv', header=TRUE, stringsAsFactors=FALSE)
C3 <- read.csv('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C3.csv', header=TRUE, stringsAsFactors=FALSE)
C4 <- read.csv('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C4.csv', header=TRUE, stringsAsFactors=FALSE)
C5 <- read.csv('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/C5.csv', header=TRUE, stringsAsFactors=FALSE)

C.tot <- rbind(C0,C1)
C.tot <- rbind(C.tot,C2)
C.tot <- rbind(C.tot,C3)
C.tot <- rbind(C.tot,C4)
C.tot <- rbind(C.tot,C5)
C.tot$ID <- as.factor(C.tot$ID)

png('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/MSMS per Duty Cycle-v2.png', width=10, height=6, units="in", res=300) 
ggplot(C.tot, aes(x=MS.scan.index, y=Freq, colour=ID))+geom_point(size=0.1)+
geom_smooth()+ xlim(0, 5000)+xlab('Duty Cycle')+ylab('Frequency')+
theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))
#+ scale_colour_discrete(breaks=c("All MSMS","Unidentified","[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
#geom_area(aes(fill=ID), alpha=0.1)+ggtitle("0")+geom_line()
dev.off()


#--------Precursor m/z Distribution

msmsSc2 <- msmsSc
msmsSc2[msmsSc$ID==0,]$ID <- "Unidentified"
msmsSc2[msmsSc$ID==1,]$ID <- "[0, 0.02]"
msmsSc2[msmsSc$ID==2,]$ID <- "(0.02, 0.05]"
msmsSc2[msmsSc$ID==3,]$ID <- "(0.05, 0.1]"
msmsSc2[msmsSc$ID==4,]$ID <- "(0.1, 1]"
msmsSc2$ID <- as.factor(msmsSc2$ID)
png('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/pecursor mz distribution for all scans(all PEPs)-v3.png', width=10, height=6, units="in", res=300)  
ggplot(msmsSc2, aes(m.z, colour = ID, show.legend = NA)) +
  geom_density(alpha=0.2, size=1.5, show.legend = NA)+xlab('Precursor m/z')+
  ylab('Probability')+labs(title='Distribution of Precursor m/z')+
  theme(text = element_text(size=20)) #Text Size
  #scale_fill_discrete(labels=c("Unidentified", "[0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "(0.1, 1]"))
dev.off()

png('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/pecursor mz distribution for all scans(iden-uniden)-v2.png', width=10, height=6, units="in", res=300)  
ggplot(msmsSc, aes(m.z, fill = Identified, colour = Identified)) +
  geom_density(alpha=0.7)+xlab('Precursor m/z')+
  ylab('Probability')+labs(title='Distribution of Precursor m/z')+
  theme(text = element_text(size=20)) #Text Size
dev.off()


#------------Computing the fractions:

C0$all.count <- C5$Freq[match(C0$MS.scan.index, C5$MS.scan.index)]
C1$all.count <- C5$Freq[match(C1$MS.scan.index, C5$MS.scan.index)]
C2$all.count <- C5$Freq[match(C2$MS.scan.index, C5$MS.scan.index)]
C3$all.count <- C5$Freq[match(C3$MS.scan.index, C5$MS.scan.index)]
C4$all.count <- C5$Freq[match(C4$MS.scan.index, C5$MS.scan.index)]

C0$Fraction <- (C0$Freq/C0$all.count)
C1$Fraction <- (C1$Freq/C1$all.count)
C2$Fraction <- (C2$Freq/C2$all.count)
C3$Fraction <- (C3$Freq/C3$all.count)
C4$Fraction <- (C4$Freq/C4$all.count)

C.tot <- rbind(C0,C1)
C.tot <- rbind(C.tot,C2)
C.tot <- rbind(C.tot,C3)
C.tot <- rbind(C.tot,C4)
C.tot$ID <- as.factor(C.tot$ID)

png('C:/Users/nslavov/Desktop/R-mz dist+msms-psm per duty cycle/MSMS per Duty Cycle-Fractions of All Scans-v2.png', width=10, height=6, units="in", res=300) 
ggplot(C.tot, aes(x=MS.scan.index, y=Fraction, colour=ID))+
  geom_point(size=0.1)+ xlim(0, 5000)+xlab('Duty Cycle')+ylab('Fractions of Scans')+
  theme(text = element_text(size=10), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20))+
  geom_smooth()
#geom_area(aes(fill=ID), alpha=0.1)+ggtitle("0")
#stat_bin(mapping = NULL, data = NULL, geom = "bar", position = "stack", width = 0.9, drop = FALSE, right = FALSE, binwidth = NULL, origin = NULL, breaks = NULL, ...)
dev.off()


#-----BINING:
#Let's say we want to create ten bins with equal number of observations in each bin
bins<-10
#The cutpoints variable holds a vector of the cutpoints used to bin the data. Finally we perform the binning itself to form the discretized variable
cutpoints<-quantile(x,(0:bins)/bins)
binned <-cut(x,cutpoints,include.lowest=TRUE)
#And there you have it. The binned vector holds our new categorical variable which can then be used for further analysis.















