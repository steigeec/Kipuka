# Set up my working environment
setwd("H:/My Drive/Kipuka/Data")
library(reshape2)
library(tidyverse)
library(vegan)
library(dplyr)

otu <- read.csv("OTUs.csv")
                   
OTU <- otu[14:nrow(otu), 29:ncol(otu)] 
rownames(OTU) <- OTU[,1]
OTU <- as.data.frame(t(OTU[,-1]))
names(OTU)[1:2] <- c("OTU", "zOTU")
# Exclude blank (NA or empty string) values 
OTUtoKeep <- OTU[grepl("OTU", OTU[[1]]), ]
OTUtoKeep[3:ncol(OTUtoKeep)] <- lapply(OTUtoKeep[3:ncol(OTUtoKeep)], as.numeric)



# OTU
OTU3 <- OTUtoKeep %>%
  group_by(OTU) %>%
  summarise(across(2:(ncol(OTUtoKeep)-1), sum, na.rm = TRUE)) # Sum columns 3 to last
OTU3 <- OTU3[,-1]
OTU3 <- as.data.frame(t(OTU3))
OTU3[] <- lapply(OTU3, as.numeric)

zOTUbeta <- vegdist(OTU3, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-as.data.frame(zOTUbeta)
dist_long <- melt(as.matrix(zOTUbeta))
write.csv(dist_long, "OTU3_Bray.csv", quote=F, row.names=F)

zOTUbeta <- vegdist(OTU3, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-as.data.frame(zOTUbeta)
dist_long <- melt(as.matrix(zOTUbeta))
write.csv(dist_long, "OTU3_jaccard.csv", quote=F, row.names=F)




# zOTU
# Identify the values in column 1 that occur more than once
values_to_keep <- names(which(table(OTUtoKeep[[1]]) > 1))

# Subset the dataframe to keep only rows where column 1 matches those values
OTUtoKeep_filtered <- OTUtoKeep[OTUtoKeep[[1]] %in% values_to_keep, ]
OTUtoKeep_filtered <- OTUtoKeep_filtered[,-c(1,2)]
OTUtoKeep_filtered <- as.data.frame(t(OTUtoKeep_filtered))
OTUtoKeep_filtered[] <- lapply(OTUtoKeep_filtered, as.numeric)

zOTUbeta <- vegdist(OTUtoKeep_filtered, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-as.data.frame(zOTUbeta)
dist_long <- melt(as.matrix(zOTUbeta))
write.csv(dist_long, "zOTU_Bray.csv", quote=F, row.names=F)

zOTUbeta <- vegdist(OTUtoKeep_filtered, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-as.data.frame(zOTUbeta)
dist_long <- melt(as.matrix(zOTUbeta))
write.csv(dist_long, "zOTU_jaccard.csv", quote=F, row.names=F)
