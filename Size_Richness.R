#Project: Kipukas               #
#Script: Size_Richness.R     #
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
#Kipuka size versus SR/SROTU

richness_mod_1 <- richness[richness$Site=="Center" | richness$Site=="Edge",]

a <- ggplot() + 
  geom_smooth(method='lm', data=richness_mod_1, aes(x=Arealog, y=SROTU, colour=Site, fill=Site), size=1, alpha=0.20)+ #linetype=variable, 
  geom_point(data=richness_mod_1,aes(x=Arealog, y=SROTU, colour=Site), alpha=0.70, size=6, stroke = 3) + #, shape=variable
  #scale_shape_manual("Site", values=c("SR" = 0, "SROTU"=15)) +
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual(values=SiteColors)+  
  facet_wrap(~Site)+                 
  #scale_linetype_discrete(values=c(2,5)) +
  labs(title="A.   Size by 3% OTU richness", x="Log area ("~m^2~")", y="3% OTU richness") +
  KipukaTheme +
  guides(color="none", fill="none") +#shape="none", 
  theme(strip.text = element_text(size = 30), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=50), 
        axis.text = element_text(size=45), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top")

b <- ggplot() + 
  geom_smooth(method='lm', data=richness_mod_1, aes(x=Arealog, y=SR, colour=Site, fill=Site), size=1, alpha=0.20)+ #linetype=variable, 
  geom_point(data=richness_mod_1,aes(x=Arealog, y=SR, colour=Site), alpha=0.70, size=6, stroke = 3) + #, shape=variable
  #scale_shape_manual("Site", values=c("SR" = 0, "SROTU"=15)) +
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual(values=SiteColors)+  
  facet_wrap(~Site)+                              
  #scale_linetype_discrete(values=c(2,5)) +
  labs(title="B.   Size by zOTU richness", x="Log area ("~m^2~")", y="zOTU richness") +
  KipukaTheme +
  guides(color="none", fill="none") +#shape="none", 
  theme(strip.text = element_text(size = 30), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=50), 
        axis.text = element_text(size=45), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top")    
                   
################################################################################
#haplotype richness within OTUs for Kipuka Centers and Kipuka edges

c <- ggplot() + 
  geom_smooth(method='lm', data=richness_mod_1,aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site, fill=Site), size=1, alpha=0.20)+
  geom_point(data=richness_mod_1,aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site), shape=18, alpha=0.70, size=6, stroke = 3) + 
  facet_wrap(~Site)+
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual(values=SiteColors, limits = c("Center", "Edge"))+                 
  labs(title="C.   Size vs haplotype richness", x="Log area ("~m^2~")", y="Haplotype richness within OTUs") +
  KipukaTheme +
  guides(color="none", fill="none") +
  theme(strip.text = element_text(size = 30), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=50), 
        axis.text = element_text(size=45), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top")

###########################################################
#Plot these three together

jpeg("Figures/Figure3.jpg", width=3000, height=1000)
plot_grid(a, b, c, ncol = 3, rel_widths = c(1, 1, 1))
dev.off()
                   