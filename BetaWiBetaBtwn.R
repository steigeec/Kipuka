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

#center versus all other types? (*1K*C)
co<-as.numeric(c())
cc<-as.numeric(c())
#edge versus all other types? (*1K*E)
eo<-as.numeric(c())
ee<-as.numeric(c())
#Stainbeck versus all other types? (HEH*K1)
so<-as.numeric(c())
ss<-as.numeric(c())
#Lava versus all other types? (*5KL*)
lo<-as.numeric(c())
ll<-as.numeric(c())
#Kona versus all other types? (*HEK*K1)
ko<-as.numeric(c())
kk<-as.numeric(c())

for (ROW in 1:nrow(zbeta)){
  site1<-zbeta$site1[ROW]
  site2<-zbeta$site2[ROW]
  if (grepl("*1K*C",site1)==TRUE){
    print("yes")
    }
  }










#Compare beta between and within site types for 3% OTU
BC3<-read.csv("3OTU.csv")
rownames(BC3)<-BC3[,1]
BC3<-BC3[,-1]
