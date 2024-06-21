
#Set up my working environment
setwd("G:/My Drive/Kipuka")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
library(vegan)
font_import()


richness <- read.csv("merged_by_site_2.csv")
dist_beta <- read.csv("Distance_v_beta.csv")
dist_diff <- read.csv("Distance_v_differentiation.csv")
geo_dist<-read.csv("geo_dist.csv")
OTU <- read.csv("OTUs.csv")

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainbeck"="#999933")

KipukaTheme <- theme(axis.title=element_text(size=30), 
        axis.text = element_text(size=25, angle=45), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        plot.title=element_text(size=30), 
        legend.text=element_text(size=25), 
        legend.key.height = unit(1, "cm"), 
        legend.key.width = unit(1.5,"cm"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.title = element_text(size=25), 
        text = element_text(family = "serif"), 
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.spacing.y = unit(0.1,"cm")) 
         
#################################################################################
#Beta diversity across orders
#For everything except hymenoptera,   calculate Bray-Curtis distances between site pairs

#wrangle geographic distances
rownames(geo_dist) <- geo_dist[,1]
geo_dist <- geo_dist[,-1]
#remove x from all the column names
names(geo_dist)<-sub("X*", "", names(geo_dist))
geo_dist<-as.matrix(geo_dist)
geo_dist<-data.frame(col=colnames(geo_dist)[col(geo_dist)], row=rownames(geo_dist)[row(geo_dist)], dist=c(geo_dist))
#make an index column that reps this particular combination of sites
geo_dist$index<-paste(geo_dist$col, geo_dist$row, sep="_")
geo_dist$log_dist<-log(geo_dist$dist+0.00001)
geo_dist<-geo_dist[,c(4, 5)]
#xy <- t(combn(colnames(geo_dist), 2))
#geo_dist <- data.frame(xy, dist=geo_dist[xy])

OTU <- OTU[8:67, 29:3051] 
#Make first column the row names
OTU["13",1]<-"ZotuID"
OTU["15",1]<-"OtuID"
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)

#What are the orders for which we need to tailor?
unique(OTU$Order)
#Classes: Acari
#Order: (1) "Coleoptera"       (2) "Diptera"          (3) "Hemiptera"        (4) "Lepidoptera"        (5) "Psocoptera"      (6) Aranea 

orders<-c("Coleoptera", "Diptera", "Hemiptera", "Lepidoptera", "Psocoptera", "Araneae")
#
for (ORDER in 1:length(orders)){
        O<-orders[ORDER]
        acari <- OTU[OTU$Order==O,] 
        #get rid of empty rows
        acari <- acari[rowSums(is.na(acari)) != ncol(acari), ]  
        acari<-acari[,9:60]
        #Sites must be rows, and species are columns
        acari<-as.data.frame(t(as.matrix(acari)))
        acari[] <- lapply(acari, as.numeric)
        acari_beta <- vegdist(acari, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
        #Convert distance matrix into longform
        acari_beta<-as.matrix(acari_beta)
        acari_beta<-data.frame(col=colnames(acari_beta)[col(acari_beta)], row=rownames(acari_beta)[row(acari_beta)], dist=c(acari_beta))
        #Add attributes of each site
        #First we add whether it's center, edge, etc etc etc
        acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="col", by.y="ï..ID")
        acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="row", by.y="ï..ID")
        #If Site.x and Site.y are not the same (e.g. not both Center and Center), then throw out that row
        acari_beta<-acari_beta[acari_beta$Site.x==acari_beta$Site.y,]
        #Now add the distances between these sites
        acari_beta$index<-paste(acari_beta$row, acari_beta$col, sep="_")
        acari_beta<-merge(acari_beta, geo_dist, by.x="index", by.y="index")
        #remove the same-site pairs
        acari_beta<-acari_beta[acari_beta$row!=acari_beta$col,]
        acari_beta$order <- O
        assign(paste0(O, "_beta"), acari_beta)  
}

# acari!!
acari <- OTU[OTU$Class=="Acari",] 
acari <- acari[rowSums(is.na(acari)) != ncol(acari), ]  
acari<-acari[,9:60]
#Sites must be rows, and species are columns
acari<-as.data.frame(t(as.matrix(acari)))
acari[] <- lapply(acari, as.numeric)
acari_beta <- vegdist(acari, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
#Convert distance matrix into longform
acari_beta<-as.matrix(acari_beta)
acari_beta<-data.frame(col=colnames(acari_beta)[col(acari_beta)], row=rownames(acari_beta)[row(acari_beta)], dist=c(acari_beta))
#Add attributes of each site
#First we add whether it's center, edge, etc etc etc
acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="col", by.y="ï..ID")
acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="row", by.y="ï..ID")
#If Site.x and Site.y are not the same (e.g. not both Center and Center), then throw out that row
acari_beta<-acari_beta[acari_beta$Site.x==acari_beta$Site.y,]
#Now add the distances between these sites
acari_beta$index<-paste(acari_beta$row, acari_beta$col, sep="_")
acari_beta<-merge(acari_beta, geo_dist, by.x="index", by.y="index")
#remove the same-site pairs
acari_beta<-acari_beta[acari_beta$row!=acari_beta$col,]
acari_beta$order <- "Acari"

#Paste all these various dataframes together
order_all<-rbind(acari_beta, Araneae_beta, Coleoptera_beta, Diptera_beta, Hemiptera_beta, Lepidoptera_beta, Psocoptera_beta) 


jpeg("Figures/Order_beta_diversity.jpg", width=1500, height=2000)
ggplot(data=order_all) + 
  geom_smooth(method='lm', aes(x=log_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20)+
  geom_point(aes(x=log_dist, y=dist, colour=Site.x), alpha=0.70, size=4, shape=18) + 
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual("Site type", values=SiteColors) +
  labs(title="Distance vs OTU beta diversity", x="Log distance (km)", y="OTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=4)+
  KipukaTheme +
  coord_cartesian(ylim=c(0.4, 1))+
  guides(colour="none")+
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=45), 
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()

#Save another version of this jpeg, now only highlighting edge and center turnover
jpeg("Figures/Order_beta_diversity_kipukas.jpg", width=1500, height=2000)
ggplot(data=order_all[order_all$Site.x==c("Center", "Edge"),]) + 
  geom_smooth(method='lm', aes(x=log_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20)+
  geom_point(aes(x=log_dist, y=dist, colour=Site.x), alpha=0.70, size=4, shape=18) + 
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual("Position in kipuka", values=SiteColors, limits = c("Center", "Edge")) +
  labs(title="Distance vs OTU beta diversity", x="Log distance (km)", y="OTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=4)+
  KipukaTheme +
  coord_cartesian(ylim=c(0.4, 1))+
  guides(colour="none")+
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=45), 
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()
#################################################################################
#Repeat the same as above, now for method jaccard

orders<-c("Coleoptera", "Diptera", "Hemiptera", "Lepidoptera", "Psocoptera", "Araneae")
#
for (ORDER in 1:length(orders)){
        O<-orders[ORDER]
        acari <- OTU[OTU$Order==O,] 
        #get rid of empty rows
        acari <- acari[rowSums(is.na(acari)) != ncol(acari), ]  
        acari<-acari[,9:60]
        #Sites must be rows, and species are columns
        acari<-as.data.frame(t(as.matrix(acari)))
        acari[] <- lapply(acari, as.numeric)
        acari_beta <- vegdist(acari, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
        #Convert distance matrix into longform
        acari_beta<-as.matrix(acari_beta)
        acari_beta<-data.frame(col=colnames(acari_beta)[col(acari_beta)], row=rownames(acari_beta)[row(acari_beta)], dist=c(acari_beta))
        #Add attributes of each site
        #First we add whether it's center, edge, etc etc etc
        acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="col", by.y="ï..ID")
        acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="row", by.y="ï..ID")
        #If Site.x and Site.y are not the same (e.g. not both Center and Center), then throw out that row
        acari_beta<-acari_beta[acari_beta$Site.x==acari_beta$Site.y,]
        #Now add the distances between these sites
        acari_beta$index<-paste(acari_beta$row, acari_beta$col, sep="_")
        acari_beta<-merge(acari_beta, geo_dist, by.x="index", by.y="index")
        #remove the same-site pairs
        acari_beta<-acari_beta[acari_beta$row!=acari_beta$col,]
        acari_beta$order <- O
        assign(paste0(O, "_beta"), acari_beta)  
}

# acari!!
acari <- OTU[OTU$Class=="Acari",] 
acari <- acari[rowSums(is.na(acari)) != ncol(acari), ]  
acari<-acari[,9:60]
#Sites must be rows, and species are columns
acari<-as.data.frame(t(as.matrix(acari)))
acari[] <- lapply(acari, as.numeric)
acari_beta <- vegdist(acari, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
#Convert distance matrix into longform
acari_beta<-as.matrix(acari_beta)
acari_beta<-data.frame(col=colnames(acari_beta)[col(acari_beta)], row=rownames(acari_beta)[row(acari_beta)], dist=c(acari_beta))
#Add attributes of each site
#First we add whether it's center, edge, etc etc etc
acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="col", by.y="ï..ID")
acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="row", by.y="ï..ID")
#If Site.x and Site.y are not the same (e.g. not both Center and Center), then throw out that row
acari_beta<-acari_beta[acari_beta$Site.x==acari_beta$Site.y,]
#Now add the distances between these sites
acari_beta$index<-paste(acari_beta$row, acari_beta$col, sep="_")
acari_beta<-merge(acari_beta, geo_dist, by.x="index", by.y="index")
#remove the same-site pairs
acari_beta<-acari_beta[acari_beta$row!=acari_beta$col,]
acari_beta$order <- "Acari"

#Paste all these various dataframes together
order_all<-rbind(acari_beta, Araneae_beta, Coleoptera_beta, Diptera_beta, Hemiptera_beta, Lepidoptera_beta, Psocoptera_beta) 


jpeg("Figures/Order_beta_diversity_JACCARD.jpg", width=1500, height=2000)
ggplot(data=order_all) + 
  geom_smooth(method='lm', aes(x=log_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20)+
  geom_point(aes(x=log_dist, y=dist, colour=Site.x), alpha=0.70, size=4, shape=18) + 
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual("Site type", values=SiteColors) +
  labs(title="Distance vs OTU beta diversity (Jaccard)", x="Log distance (km)", y="OTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=4)+
  KipukaTheme +
  coord_cartesian(ylim=c(0.6, 1))+
  guides(colour="none")+
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=45), 
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()


#Save an alternate version of this jpeg, highlighting only kipuka data
jpeg("Figures/Order_beta_diversity_JACCARD_kipukas.jpg", width=1500, height=2000)
ggplot(data=order_all[order_all$Site.x==c("Center", "Edge"),]) + 
  geom_smooth(method='lm', aes(x=log_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20)+
  geom_point(aes(x=log_dist, y=dist, colour=Site.x), alpha=0.70, size=4, shape=18) + 
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual("Position in kipuka", values=SiteColors, limits = c("Center", "Edge")) +
  labs(title="Distance vs OTU beta diversity (Jaccard)", x="Log distance (km)", y="OTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=4)+
  KipukaTheme +
  coord_cartesian(ylim=c(0.6, 1))+
  guides(colour="none")+
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=45), 
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()
