#Project: Kipukas               #
#Script: PlottingRichness.R     #
#Author: Emma Steigerwald       #
#Date:17 Feb 2022               #
#################################

#Set up my working environment
setwd("G:/My Drive/Kipuka")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
library(vegan)
library(scales)
font_import()



#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainback"="#999933")
#Establish some themes up top to apply to all
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

#0% beta diversity
dist_beta_0 <- read.csv("Distance_v_beta_1.csv")
#zOTU beta diversity
otu <- read.csv("Distances_Without_Kona8andsmall_kipuka.csv")
geo_dist<-read.csv("geo_dist.csv")
richness <- read.csv("merged_by_site_2.csv")
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))


#################################################################################
#Beta diversity vesus distance

#Wrange geo distances
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

#Grab zOTU beta diversity metrics
acari_beta<-as.matrix(otu)
#Remove first 4 cols (unclear what they are)
acari_beta<-acari_beta[,5:10]
#Remove entirely NA rows
acari_beta <- acari_beta[rowSums(is.na(acari_beta)) != ncol(acari_beta), ] 

#Col names need to be formatted so that they match color scheme
colnames(acari_beta) <- c("Geo", "log.1", "Kona", "Stainback", "Center", "Edge")
acari_beta<-as.data.frame(acari_beta[,c(1, 3, 4, 5, 6)])

#Now rearrange from wide to long format
acari_beta <- melt(acari_beta, id= c("Geo"))

#Assign to "kipuka" or "forest"
acari_beta$group[acari_beta$variable=="Center" | acari_beta$variable=="Edge"] <- "Kipuka"
acari_beta$group[acari_beta$variable=="Kona" | acari_beta$variable=="Stainback"] <- "Continuous forest"
acari_beta<-acari_beta[complete.cases(acari_beta),]
acari_beta<-acari_beta[!duplicated(acari_beta),]

b <- ggplot(data=acari_beta) + 
  geom_point(aes(x=Geo, y=value, colour=variable), alpha=0.70, size=8, shape=19) + 
  geom_smooth(method='lm', aes(x=Geo, y=value, colour=variable), size=1, alpha=0)+ #, linetype=site
  facet_wrap(~group, ncol=2)+
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge", "Kona", "Stainback")) +
  labs(title="B.", x="Distance (km)", y="zOTU beta diversity") +
  KipukaTheme +
  guides(color = guide_legend(title = "Sites", nrow=1)) +
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +
  theme(strip.text=element_text(size=45), 
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
       legend.position = "bottom")
                           
#################################################################################
#now 3%OTU
                     
#order acari_beta by site then value... do the same with dist_beta_0                  
acari_beta$variable  <- factor(acari_beta$variable , levels=c("Center", "Edge", "Stainback", "Kona", "Lava"))                   
acari_beta <- acari_beta[with(acari_beta, order(variable, value)),] 
dist_beta_0 <- dist_beta_0[dist_beta_0$Site!="Lava",]                     
dist_beta_0$Site  <- factor(dist_beta_0$Site , levels=c("Center", "Edge", "Stainback", "Kona"))                                                          
dist_beta_0 <- dist_beta_0[with(dist_beta_0, order(Site, beta)),]                     

#join the two... 
acari_beta$threeotu <- dist_beta_0$beta                    
                     
a <- ggplot(data=acari_beta) + 
  geom_smooth(method='lm', aes(x=Geo, y=threeotu, colour=variable, fill=variable, linetype=variable), size=1, alpha=0)+
  geom_point(aes(x=Geo, y=threeotu, colour=variable), alpha=0.70, size=8, shape=19) + 
  facet_wrap(~group, ncol=2)+
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual(values=SiteColors) +
  labs(title="A.", x="Distance (km)", y="3% OTU beta diversity") +
  KipukaTheme +
  guides(color="none", shape="none", fill ="none", linetype="none") +
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +
  theme(strip.text=element_text(size=45),
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
       legend.position = "bottom")
                                
                     
#################################################################################                         
jpeg("Figures/Fig_4.jpg", width=2000, height=2000)   
plot_grid(a, b, nrow = 2, rel_heights = c(1, 1.3))                         
dev.off()
                     
