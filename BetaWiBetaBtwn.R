# BetaWiBetaBtwn.R    #
# Emma Steigerwald    #
# 17 March 2023       #
#######################

#Set up my working environment
setwd("G:/My Drive/Kipuka")
library(extrafont)
library(reshape2)
library(tidyverse)
library(scales)

#Compare beta between and within site types for zOTU
BC<-read.csv("zotu.csv")
rownames(BC)<-BC[,1]
BC<-BC[,-1]
zbeta<-as.data.frame(matrix(0, nrow = nrow(BC)*nrow(BC), ncol = 3))
names(zbeta)<-c("site1", "site2", "beta")
for (ROW in 1:nrow(BC)){
  Rname<-row.names(BC)[ROW]
  for (COL in 1:ncol(BC)){
    Cname<-names(BC)[COL]
    distance<-BC[ROW,COL]
    #The row # to place in the output matrix will be row# + (col#-1)*nrow(BC)
    zbeta$beta[ROW+((COL-1)*nrow(BC))]<-distance
    zbeta$site1[ROW+((COL-1)*nrow(BC))]<-Rname
    zbeta$site2[ROW+((COL-1)*nrow(BC))]<-Cname
  }
}

#











#Compare beta between and within site types for 3% OTU
BC3<-read.csv("3OTU.csv")
rownames(BC3)<-BC3[,1]
BC3<-BC3[,-1]
