Set up my working environment
setwd("H:/My Drive/Kipuka/Data")
library(reshape2)
library(tidyverse)
library(vegan)

otu <- read.csv("OTUs.csv")
                   
OTU <- otu[16:nrow(otu), 29:ncol(otu)] 
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)
#get rid of empty rows
#OTU <- OTU[rowSums(is.na(OTU)) != ncol(OTU), ]  
OTU<-as.data.frame(t(as.matrix(OTU)))
OTU[] <- lapply(OTU, as.numeric)

zOTUbeta <- vegdist(OTU, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-data.frame(col=colnames(zOTUbeta)[col(zOTUbeta)], row=rownames(zOTUbeta)[row(zOTUbeta)], dist=c(zOTUbeta))

write.csv(zOTUbeta, "zOTU_Bray.csv", quote=F, row.names=F)
