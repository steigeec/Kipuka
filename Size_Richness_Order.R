#Project: Kipukas                  #
#Script: Size_Richness_Order.R     #
#Author: Emma Steigerwald          #
#Date:17 Feb 2022                  #
####################################

#Set up my working environment
setwd("G:/My Drive/Kipuka")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
library(vegan)
font_import()
library(scales)


richness <- read.csv("merged_by_site_2.csv")
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))
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
#Let's just try all the orders! 

richness_mod_0 <- richness[richness$Site=="Center" | richness$Site=="Edge",]
richness_mod_0 <- melt(richness, idvars = c("SiteID", "Area"), measure = c("Araneae", "Pscoptera", "Hemiptera", "Lepidoptera", "Coleoptera", "Diptera"))
richness_mod_0$Area <- round(richness_mod_0$Area, 10)

#Put in correct order
#richness_mod_0$Area <- factor(richness_mod_0$Area, levels=c("Lava", "3", "4", "5", "Stainbeck", "Kona"))

jpeg("Figures/Order_Richness_continuous.jpg", width=3000, height=2000)
ggplot() + 
  geom_point(data=richness_mod_0,aes(x=Area, y=value, colour=Site, fill=Site), size=4.5, shape=0, stroke=2)+
  geom_smooth(method='lm', data=richness_mod_0, aes(x=Area, y=value, colour=Site, fill=Site), size=1, alpha=0.20)+
  scale_fill_manual(values=SiteColors, limits=c("Center", "Edge")) +
  scale_colour_manual(values=SiteColors) +
  facet_wrap(~variable, nrow=2, scales="free") +
  guides(fill=guide_legend(nrow=1), colour="none") +
  labs(title="", x="Kipuka area ("~m^2~")", y="zOTU richness") +
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  + 
  KipukaTheme +
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title.x=element_text(size=55, margin=margin(-15,0,0,0)), 
       axis.title.y=element_text(size=55, margin=margin(0,-15,0,0)), 
        axis.text.y = element_text(size=50, angle=45), 
        axis.text.x = element_text(size=50, angle=45, vjust=0.5, hjust=0.5), 
        plot.title=element_text(size=55), 
        legend.text=element_text(size=50), 
        legend.title = element_text(size=50),
       legend.position = "top", 
        plot.margin = margin(1,1,.01,1, "cm"))
dev.off()


                         











#############################OLD VERSION OF PLOTS################################################################
#################################################################################
#How is Araenea (predator!!!) richness impacted by fragment size as opposed to Pscoptera (barklice-- scavenger/detritovore) 
#richness_mod_0 <- richness[richness$Site=="Center" | richness$Site=="Edge",]
richness_mod_0 <- melt(richness, idvars = c("SiteID", "Arealog"), measure = c("Araneae", "Pscoptera"))
richness_mod_0$Arealog <- round(richness_mod_0$Arealog, 0)
richness_mod_0$Arealog[richness_mod_0$Site=="Lava" & is.na(richness_mod_0$Arealog)] <- "Lava"
richness_mod_0$Arealog[richness_mod_0$Site=="Kona" & is.na(richness_mod_0$Arealog)] <- "Kona"
richness_mod_0$Arealog[richness_mod_0$Site=="Stainbeck" & is.na(richness_mod_0$Arealog)] <- "Stainbeck"

#Put in correct order
richness_mod_0$Arealog <- factor(richness_mod_0$Arealog, levels=c("Lava", "3", "4", "5", "Stainbeck", "Kona"))

jpeg("Figures/PredatorVsDetritovore.jpg", width=1000, height=1000)
ggplot() + 
  geom_boxplot(data=richness_mod_0,aes(x=Arealog, y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  facet_wrap(~variable, nrow=1, scales="free") +
  guides(fill=guide_legend(nrow=2)) +
  labs(title="Predator v scavenger richness", x="Log area ("~km^2~")", y="Species richness") +
  KipukaTheme +
  theme(strip.text = element_text(size = 30), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=45), 
        axis.text = element_text(size=40, angle=45), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()

#################################################################################
#Let's just try all the orders! 
#richness_mod_0 <- richness[richness$Site=="Center" | richness$Site=="Edge",]
richness_mod_0 <- melt(richness, idvars = c("SiteID", "Arealog"), measure = c("Araneae", "Pscoptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Acari", "Coleoptera", "Diptera"))
richness_mod_0$Arealog <- round(richness_mod_0$Arealog, 0)
richness_mod_0$Arealog[richness_mod_0$Site=="Lava" & is.na(richness_mod_0$Arealog)] <- "Lava"
richness_mod_0$Arealog[richness_mod_0$Site=="Kona" & is.na(richness_mod_0$Arealog)] <- "Kona"
richness_mod_0$Arealog[richness_mod_0$Site=="Stainbeck" & is.na(richness_mod_0$Arealog)] <- "Stainbeck"

#Put in correct order
richness_mod_0$Arealog <- factor(richness_mod_0$Arealog, levels=c("Lava", "3", "4", "5", "Stainbeck", "Kona"))

jpeg("Figures/Order_Richness_1.jpg", width=4000, height=2000)
ggplot() + 
  geom_boxplot(data=richness_mod_0,aes(x=Arealog, y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  facet_wrap(~variable, nrow=2, scales="free") +
  guides(fill=guide_legend(nrow=2)) +
  labs(title="Predator v scavenger richness", x="Log area ("~km^2~")", y="Species richness") +
  KipukaTheme +
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=55), 
        axis.text.y = element_text(size=50, angle=45), 
        axis.text.x = element_text(size=50, angle=45, vjust=0), 
        plot.title=element_text(size=55), 
        legend.text=element_text(size=50), 
        legend.title = element_text(size=50),
       legend.position = "top", 
        plot.margin = margin(1,1,.01,1, "cm"))
dev.off()

