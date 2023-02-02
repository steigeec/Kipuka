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
library(scales)


richness <- read.csv("merged_by_site_2.csv")
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainbeck"="#999933")
#Establish some themes up top to apply to all
KipukaTheme <- theme(axis.title=element_text(size=30), 
        axis.text = element_text(size=25, angle=45), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        plot.title=element_text(size=30), 
        legend.text=element_text(size=25), 
        legend.spacing.x = unit(1.0, 'cm'),
        legend.key.height = unit(1, "cm"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.title = element_text(size=25), 
        text = element_text(face="plain", family="Calibri"),
        legend.box.background = element_rect(fill = "white", color = "black")) 
                   
#################################################################################
#Kipuka size versus SR/SROTU
richness_mod_2 <- melt(richness, idvars = c("SiteID", "Site"), measure = c("SR", "SROTU"))
richness_mod_2 <- richness_mod_2[order(richness_mod_2$value, decreasing = TRUE),]  

# New facet label names for supp variable
supp.labs <- c("zOTU richness", "3% OTU richness")
names(supp.labs) <- c("SR", "SROTU")                   

#Reorder facets
richness_mod_2$variable <- factor(richness_mod_2$variable, levels = rev(c("SR", "SROTU")))                     
                               
a <- ggplot() + 
  geom_boxplot(data=richness_mod_2,aes(x=reorder(Site, value), y=value, fill=Site), color="black", size=1)+
  facet_wrap(~variable, scales="free", 
             labeller = labeller(variable = supp.labs)) +
  scale_fill_manual(values=SiteColors) +
  labs(title="A.", x="") +
  KipukaTheme +
  theme(strip.text = element_text(size = 35), 
        axis.text.y = element_text(angle=45, size=45), 
        axis.text.x = element_text(angle=45, size=45, vjust=0.6), 
        axis.title.y = element_blank(), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=50), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=45, hjust=0.4), 
        legend.title = element_blank(),
        legend.key.width = unit(4,"cm"), 
       legend.position = "top", 
        plot.margin = margin(0.1,-1.5,2.5,2, "cm"))+ guides(shape = guide_legend(override.aes = list(size = 4)))

CenterzOTU <- lm(richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SR"]~richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SR"])
Center3otu<-lm(richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SROTU"]~richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SROTU"])
EdgezOTU<-lm(richness_mod_2$value[richness_mod_2$variable=="SR" & richness_mod_2$Site=="Edge"]~richness_mod_2$value[richness_mod_2$variable=="SR" & richness_mod_2$Site=="Edge"])
Edge3otu<-lm(richness_mod_2$value[richness_mod_2$variable=="SROTU" & richness_mod_2$Site=="Edge"]~richness_mod_2$value[richness_mod_2$variable=="SROTU" & richness_mod_2$Site=="Edge"])

m <- lm(y ~ x)
r2 = format(summary(m)$r.squared, digits = 3)))                     
                  

b<-ggplot() + 
  geom_smooth(method='lm', data=richness_mod_2[richness_mod_2$Site=="Center" | richness_mod_2$Site=="Edge",], aes(x=Area, y=value, colour=Site, fill=Site, linetype=variable), size=1, alpha=0.20)+  
  geom_point(data=richness_mod_2[richness_mod_2$Site=="Center" | richness_mod_2$Site=="Edge",],aes(x=Area, y=value, colour=Site, shape=variable), alpha=0.70, size=6, stroke = 3) + 
  scale_shape_manual("Site", values=c("SR" = 0, "SROTU"=15), labels=c("SR"="zOTU","SROTU"="3% OTU")) +
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual(values=SiteColors)+ 
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +                 
  facet_wrap(~Site)+                 
  labs(title="B.", x="Kipuka area ("~m^2~")", y="OTU richness") +
  KipukaTheme +
  guides(color="none", fill="none", linetype="none") + 
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title.y=element_text(size=50, vjust=-0.5, hjust=0.5), 
       axis.title.x=element_text(size=50, vjust=2, hjust=0.5), 
        axis.text.y = element_text(size=45), 
        axis.text.x = element_text(size=45, vjust=1), 
        plot.title=element_text(size=50), 
        legend.key.width = unit(7,"cm"), 
        legend.text=element_text(size=45, hjust=0.5), 
        legend.title = element_blank(),
       legend.position = "top", 
        plot.margin = margin(0.1,2.5,2.5,3, "cm"))+ guides(shape = guide_legend(override.aes = list(size = 4)))

jpeg("Figures/Figure3.jpg", width=2000, height=1000)
plot_grid(a, b, ncol = 2, rel_widths = c(1, 2))
dev.off()                     
 
                     

                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
##########################################
#old plot format                     
b <- ggplot() + 
  geom_smooth(method='lm', data=richness[richness$Site=="Center" | richness$Site=="Edge",], aes(x=Arealog, y=SR, colour=Site, fill=Site), size=1, alpha=0.20)+ #linetype=variable, 
  geom_point(data=richness[richness$Site=="Center" | richness$Site=="Edge",],aes(x=Arealog, y=SR, colour=Site), alpha=0.70, size=6, stroke = 3) + #, shape=variable
  geom_hline(yintercept=mean(richness$SR[richness$Site=="Kona"]), colour="#117733", lwd=3)+
  geom_hline(yintercept=mean(richness$SR[richness$Site=="Stainbeck"]), colour="#999933", lwd=3)+
  geom_hline(yintercept=mean(richness$SR[richness$Site=="Lava"]), colour="#888888", lwd=3, alpha=0.6)+
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge", "Kona", "Stainbeck")) +
  scale_fill_manual(values=SiteColors)+  
  facet_wrap(~Site)+                              
  #scale_linetype_discrete(values=c(2,5)) +
  labs(title="B.", x="Log area ("~m^2~")", y="zOTU richness") +
  KipukaTheme +
  guides(color="none", fill="none") +#shape="none", 
  theme(strip.text = element_text(size = 35), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=50), 
        axis.text.y = element_text(size=45), 
        axis.text.x = element_text(size=45, vjust=0.6),  
        plot.title=element_text(size=50), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top")    
                   
################################################################################
#haplotype richness within OTUs for Kipuka Centers and Kipuka edges

c <- ggplot() + 
  geom_smooth(method='lm', data=richness[richness$Site=="Center" | richness$Site=="Edge",], aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site, fill=Site), size=1, alpha=0.20)+ #linetype=variable, 
  geom_point(data=richness[richness$Site=="Center" | richness$Site=="Edge",],aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site), alpha=0.70, size=6, stroke = 3) + #, shape=variable
  geom_hline(yintercept=mean(richness$HaplotypeRichnessWithin[richness$Site=="Kona"]), colour="#117733", lwd=3)+
  geom_hline(yintercept=mean(richness$HaplotypeRichnessWithin[richness$Site=="Stainbeck"]), colour="#999933", lwd=3)+
  geom_hline(yintercept=mean(richness$HaplotypeRichnessWithin[richness$Site=="Lava"]), colour="#888888", lwd=3, alpha=0.6)+
  facet_wrap(~Site)+
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge", "Kona", "Stainbeck")) +
  scale_fill_manual(values=SiteColors, limits = c("Center", "Edge"))+                 
  labs(title="C.", x="Log area ("~m^2~")", y="Haplotype richness within OTUs") +
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

#jpeg("Figures/Figure3.jpg", width=2000, height=1000)
#plot_grid(a, b, ncol = 2, rel_widths = c(1, 1))
#dev.off()
                   
