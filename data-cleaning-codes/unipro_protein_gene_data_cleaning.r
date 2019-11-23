pg <- read.delim('./uniprot_Gene-Prot.txt', header=TRUE, sep = "\t")
pg <- data.frame(lapply(pg, as.character), stringsAsFactors=FALSE) #strsplit works only with characters

#Assuming the first element is always the 'key' and the rest are aliases, split the gene names,
#identify the key, then group all aliases by key, and standardize each element to include aliases:
elts <- strsplit(pg$Gene.names, " ")
keys <- sapply(elts, "[[", 1)
values <- split(unlist(elts), rep(keys, lengths(elts)))
pg$Gene.names <- lapply(values, unique)[keys]

#Use the length of each standardized gene name to replicate the entry names,
#and match these to the unlisted, split gene names:
pgsplt <- data.frame(Entry.name = rep(pg$Entry.name, lengths(pg$Gene.names)), Gen.name = unlist(pg$Gene.names))

#Getting the frequency of genes (Which shows the number of associated proteins):
pgfrq <- data.frame (table(pgsplt$Gen.name))

#Adding the count of proteins of each gene in a column at the end of UTR Gene list:
C$prot.count <- pgfrq$Freq[match(C$NULL.Gene, pgfrq$Var1)]

write.csv(C, './Dropbox/Genetics/R/New Results/5UTR-Peak-Comp-w prot.count.csv')

write.csv(pgfrq, './Dropbox/Genetics/R/New Results/p-g-freq.csv')

Z <- data.frame(V1 = c("A", "B", "C", "D", "E"), V2 = c('1', '3', '5', '2', '4'))
Zt <- data.frame(V1 = c("A", "A", "B", "D", "E"))
Zt$V2 <- Z$V2[match(Zt$V1, Z$V1)]
