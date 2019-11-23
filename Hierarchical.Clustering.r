prot <- read.csv('C:/Users/nslavov/Desktop/fwdpotentialcollaboration/ORGalexander_tmt10_20160217_table.csv', header = TRUE, stringsAsFactors = FALSE)
clusters <- hclust(dist(prot$PMA.1/mean(prot$PMA.1), method = "euclidean"))
PMA <- data.frame(C1=as.numeric(), C2=as.numeric())

C1 <- as.numeric(prot$PMA.1/mean(prot$PMA.1))
C2 <- as.numeric(prot$PMA.2/mean(prot$PMA.2))
PMA <- data.frame(C1, C2)

PMA.clusters <- hclust(dist(PMA, method = "euclidean"))

plot(PMA.clusters)

#we can cut off the tree at the desired number of clusters using cutree.

clusterCut <- cutree(PMA.clusters, 4)
plot(clusterCut)
install.packages('ape')
library('ape')
plot(as.phylo(clusterCut), type = "fan")

table(clusterCut, prot$Protein_name)

install.packages('gplots')
library('gplots')
clusters <- hclust(dist(iris[, 3:4]))
plot(clusterCut)

library('ggplot2')

ggplot(PMA, aes(C1, C2)) + 
  geom_point(alpha = 0.4, size = 3.5) + geom_point(col = clusterCut) + 
  scale_color_manual(values = c('black', 'red', 'green', 'blue'))

prot.2 <- prot[,5:25]


heatmap(r.matrix, distfun=dist, hclustfun=function(d) hclust(d, method="ward"))

heatmap(as.numeric(prot[,6:25]), distfun = dist(prot[,6:25], method = "euclidean"), hclustfun = hclust(dist(prot[,6:25], method = "euclidean")))
scaled <- scale(States[,-1]) # scale all but the first column to make information comparable
heatmap.2(scaled, # specify the (scaled) data to be used in the heatmap
          cexRow=0.5, cexCol=0.95, # decrease font size of row/column labels
          scale="none", # we have already scaled the data
          trace="none") # cleaner heatmap

for (i in 2:21) {
  prot.2[,i] <- prot.2[,i]-min(prot.2[,i])/max(prot.2[,i])-min(prot.2[,i])
}

heatmap.2(as.matrix(prot.2[,2:21]), distfun=function(d) dist(d, method="euclidean"), hclustfun=hclust)

h <- hclust(dist(prot.2[,2:21], method = "euclidean"))
plot(h)


















