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
font_import()


richness <- read.csv("merged_by_site_2.csv")

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#999933", "Stainbeck"="#117733")
#Establish some themes up top to apply to all
KipukaTheme <- theme(axis.title=element_text(size=30), 
        axis.text = element_text(size=25), 
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
        text = element_text(family = "serif")) 

#################################################################################
#Kipuka size versus SR/SROTU

richness_mod <- richness[richness$Site=="Center" | richness$Site=="Edge",]
richness_mod <- melt(richness_mod, idvars = c("SiteID", "Arealog"), measure = c("SR", "SROTU"))


a <- ggplot() + 
  geom_smooth(method='lm', data=richness_mod,aes(x=Arealog, y=value, colour=Site, linetype=variable), size=1, alpha=0.20)+
  geom_point(data=richness_mod,aes(x=Arealog, y=value, colour=Site, shape=variable), alpha=0.70, size=6, stroke = 3) + 
  scale_shape_manual("Site", values=c("SR" = 0, "SROTU"=15)) +
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  #scale_linetype_discrete(values=c(2,5)) +
  labs(title="Size vs species richness", x="Log area ("~km^2~")", y="Species richness") +
  KipukaTheme +
  theme(panel.grid.major = element_line(
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
        legend.title = element_text(size=40))

################################################################################
#haplotype richness within OTUs for Kipuka Centers and Kipuka edges

richness_mod_1 <- richness[richness$Site=="Center" | richness$Site=="Edge",]

b <- ggplot() + 
  geom_smooth(method='lm', data=richness_mod_1,aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site), size=1, alpha=0.20)+
  geom_point(data=richness_mod_1,aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site), shape=18, alpha=0.70, size=6, stroke = 3) + 
  facet_wrap(~Site)+
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  labs(title="Size vs haplotype richness within OTUs", x="Log area ("~km^2~")", y="Haplotype richness within OTUs") +
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
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40))

################################################################################
#OTu richness, zOTU richness and haplotype diversity in different site types

richness_mod_2 <- melt(richness, idvars = c("SiteID", "Site"), measure = c("SR", "SROTU", "HaplotypeRichnessWithin"))
richness_mod_2 <- richness_mod_2[order(richness_mod_2$value, decreasing = TRUE),]  

c <- ggplot() + 
  geom_boxplot(data=richness_mod_2,aes(x=reorder(Site, value), y=value, fill=Site), color="black", size=1)+
  facet_wrap(~variable, scales="free")+
  scale_fill_manual(values=SiteColors) +
  labs(title="Size vs haplotype richness within OTUs", x="") +
  KipukaTheme +
  theme(strip.text = element_text(size = 30), 
        axis.text = element_text(angle=45), 
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
        legend.title = element_text(size=40))

###########################################################
#Plot these three together

jpeg("Figures/Figure2.jpg", width=3000, height=1000)
plot_grid(a, b, c, ncol = 3, rel_widths = c(1, 1, 1.5))
dev.off()

#################################################################################
#NMDS plot

#Create an NMDS plot with columns MDS1 and MDS2
richness_mod <- richness
richness_mod$Arealog[richness_mod$Site=="Lava" & is.na(richness_mod$Arealog)] <- 3
richness_mod$Arealog[richness_mod$Site=="Kona" & is.na(richness_mod$Arealog)] <- 9
richness_mod$Arealog[richness_mod$Site=="Stainbeck" & is.na(richness_mod$Arealog)] <- 9
richness_mod$Arealog<-richness_mod$Arealog^1.5

jpeg("Figures/NMDS_1.jpg", width=1000, height=1000)
ggplot() + 
  geom_point(data=richness_mod,aes(x=MDS1,y=MDS2,colour=Site, size=(Arealog), shape=Site), alpha=0.70) + 
  #geom_polygon(data=hull.data,aes(x=MDS1,y=MDS2,fill=grp,group=grp),alpha=0.30) + # add the convex hulls
  scale_colour_manual(values=SiteColors) +
  scale_shape_manual("Site", values=c("Center" = 16, "Edge" = 16, "Lava"=3, "Kona"=16, "Stainbeck"=16)) +
  scale_size_continuous("Log area ("~km^2~")", labels = c("3", "4", "5", "6", "7", "Continuous forest"), range=c(5.20,27), breaks=c(5.20, 8, 13.13, 14.70, 18.52, 27)) +
  labs(title="NMDS plot", x="NMDS1", y="NMDS2") +
  coord_equal() +
  scale_y_continuous(limits=c(-1,1.20)) +
  guides(colour = guide_legend(override.aes = list(size=10))) + 
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

