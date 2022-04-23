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


dist_beta_0 <- read.csv("Distance_v_beta.csv")
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

#Grab 3% OTU beta diversity metrics
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


a <- ggplot(data=acari_beta) + 
  geom_point(aes(x=Geo, y=value, colour=variable), alpha=0.70, size=8, shape=19) + 
  geom_smooth(method='lm', aes(x=Geo, y=value, colour=variable), size=1, alpha=0)+ #, linetype=site
  facet_wrap(~group, ncol=2)+
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge", "Kona", "Stainbeck")) +
  labs(title="A.", x="Distance (km)", y="3% OTU beta diversity") +
  KipukaTheme +
  guides(color = guide_legend(title = "Sites")) +
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
       legend.position = "right")
                           
#################################################################################
#Beta diversity vesus distance

#Maybe join together tables of each part...?
dist_beta <- dist_beta_0[,1:6]
i<-c(3, 4, 5, 6)
dist_beta[ , i] <- apply(dist_beta[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
names(dist_beta)[4]<-"Stainback"
                         
Center <- dist_beta[!is.na(dist_beta$Center),]
Center$site <- "Center"
Center <- dplyr::rename(Center, beta = Center)
Center <- Center[ , colSums(is.na(Center)) < nrow(Center)]                    
                         
Stainback <- dist_beta[!is.na(dist_beta$Stainback),]
Stainback$site <- "Stainback"  
Stainback <- dplyr::rename(Stainback, beta = Stainback)
Stainback <- Stainback[ , colSums(is.na(Stainback)) < nrow(Stainback)]   
                         
Kona <- dist_beta[!is.na(dist_beta$Kona),]
Kona$site <- "Kona"
Kona <- dplyr::rename(Kona, beta = Kona)
Kona <- Kona[ , colSums(is.na(Kona)) < nrow(Kona)]   
                         
Edge <- dist_beta[!is.na(dist_beta$Edge),]
Edge$site <- "Edge"
Edge <- dplyr::rename(Edge, beta = Edge)
Edge <- Edge[ , colSums(is.na(Edge)) < nrow(Edge)] 
                         
dist_beta <- rbind(Center, Stainback, Kona, Edge)
dist_beta <- dplyr::rename(dist_beta, dist = ï..dist)                         

#Assign to "kipuka" or "forest"
dist_beta$group[dist_beta$site=="Center" | dist_beta$site=="Edge"] <- "Kipuka"
dist_beta$group[dist_beta$site=="Kona" | dist_beta$site=="Stainback"] <- "Continuous forest"                         

#Remove smallest kipukas
#To do so, I remove any "Kipuka" group rows that do not match geographical distance pairs in my groomed dataset (called otu here)                       
#my_set <- as.numeric(acari_beta$Geo)
#Only keep dist_beta registers is (1) they are in continuous forest or (2) they are in that set of distance pairs                         
#to_plot_a <- dist_beta[dist_beta$group=="Continuous forest",]                         
#to_plot_b <- dist_beta[dist_beta$dist %in% my_set,]  
#to_plot <- rbind(to_plot_a, to_plot_b)                         
                         
b <- ggplot(data=dist_beta) + 
  geom_smooth(method='lm', aes(x=logidst, y=beta, colour=site, fill=site, linetype=site), size=1, alpha=0)+
  geom_point(aes(x=logidst, y=beta, colour=site), alpha=0.70, size=8, shape=19) + 
  facet_wrap(~group, ncol=2)+
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual(values=SiteColors) +
  labs(title="B.", x="Distance (km)", y="zOTU beta diversity") +
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
       legend.position = "right")
                                
                     
#################################################################################                         
jpeg("Figures/Fig_4.jpg", width=2000, height=1000)   
plot_grid(a, b, nrow = 2, rel_heights = c(1, 1))                         
dev.off()
                     
