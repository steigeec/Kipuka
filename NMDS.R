#Project: Kipukas               #
#Script: NMDS.R                 #
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
font_import()


richness <- read.csv("merged_by_site_2.csv")
dist_beta <- read.csv("Distance_v_beta.csv")
dist_diff <- read.csv("Distance_v_differentiation.csv")
geo_dist<-read.csv("geo_dist.csv")
OTU <- read.csv("OTUs.csv")

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainbeck"="#999933")
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

                   
#################################################################################
#NMDS plot

#Create an NMDS plot with columns MDS1 and MDS2
richness_mod <- richness
richness_mod$Area<-as.numeric(gsub(",","",as.character(richness_mod$Area)))
richness_mod$pointsize<-round(sqrt(richness_mod$Area)/10,0)
richness_mod$pointsize[richness_mod$Site=="Lava" & is.na(richness_mod$pointsize)] <- 2
richness_mod$pointsize[richness_mod$Site=="Kona" & is.na(richness_mod$pointsize)] <- 2
richness_mod$pointsize[richness_mod$Site=="Stainbeck" & is.na(richness_mod$pointsize)] <- 2

jpeg("Figures/NMDS_4.jpg", width=1000, height=1000)
ggplot() + 
  geom_point(data=richness_mod,aes(x=MDS1,y=MDS2,colour=Site, size=pointsize, shape=Site), alpha=0.70, stroke=3) + 
  scale_colour_manual(values=SiteColors) +
  scale_shape_manual("Site", values=c("Center" = 16, "Edge" = 16, "Lava"=3, "Kona"=2, "Stainbeck"=2)) +
  scale_size_continuous("Kipuka area ("~m^2~")", range=c(2,32), breaks=seq(2,32,5), labels=round((10*seq(2,32,5))^2,100)) +
  labs(title="", x="NMDS1", y="NMDS2") +
  coord_equal() +
  scale_y_continuous(limits=c(-1,1.20)) +
  guides(colour = guide_legend(override.aes = list(size=4))) + 
  KipukaTheme +
  theme(panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5))
dev.off()

