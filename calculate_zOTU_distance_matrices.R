# Set up my working environment
setwd("H:/My Drive/Kipuka/Data")
library(reshape2)
library(tidyverse)
library(vegan)


otu <- read.csv("OTUs.csv")
                   
OTU <- otu[14:nrow(otu), 29:ncol(otu)] 
rownames(OTU) <- OTU[,1]
OTU <- as.data.frame(t(OTU[,-1]))
names(OTU)[1:2] <- c("OTU", "zOTU")

# Exclude blank (NA or empty string) values 
OTUtoKeep <- OTU[grepl("OTU", OTU[[1]]), ]
# Identify the values in column 1 that occur more than once
values_to_keep <- names(which(table(OTUtoKeep[[1]]) > 1))

# Subset the dataframe to keep only rows where column 1 matches those values
OTUtoKeep_filtered <- OTUtoKeep[OTUtoKeep[[1]] %in% values_to_keep, ]
OTUtoKeep_filtered <- OTUtoKeep_filtered[,-c(1,2)]
OTUtoKeep_filtered <- as.data.frame(lapply(OTUtoKeep_filtered, as.numeric))

zOTUbeta <- vegdist(OTUtoKeep_filtered, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-data.frame(col=colnames(zOTUbeta)[col(zOTUbeta)], row=rownames(zOTUbeta)[row(zOTUbeta)], dist=c(zOTUbeta))

write.csv(zOTUbeta, "zOTU_Bray.csv", quote=F, row.names=F)

zOTUbeta <- vegdist(OTUtoKeep_filtered, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-data.frame(col=colnames(zOTUbeta)[col(zOTUbeta)], row=rownames(zOTUbeta)[row(zOTUbeta)], dist=c(zOTUbeta))

write.csv(zOTUbeta, "zOTU_jaccard.csv", quote=F, row.names=F)
