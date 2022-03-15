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
#Beta diversity vesus distance

#Maybe join together tables of each part...?
dist_beta <- dist_beta[1:6]
i<-c(3, 4, 5, 6)
dist_beta[ , i] <- apply(dist_beta[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))

Center <- dist_beta[!is.na(dist_beta$Center),]
Center$site <- "Center"
Center <- rename(Center, beta = Center)
Center <- Center[ , colSums(is.na(Center)) < nrow(Center)]                    
                         
Stainbeck <- dist_beta[!is.na(dist_beta$Stainbeck),]
Stainbeck$site <- "Stainbeck"  
Stainbeck <- rename(Stainbeck, beta = Stainbeck)
Stainbeck <- Stainbeck[ , colSums(is.na(Stainbeck)) < nrow(Stainbeck)]   
                         
Kona <- dist_beta[!is.na(dist_beta$Kona),]
Kona$site <- "Kona"
Kona <- rename(Kona, beta = Kona)
Kona <- Kona[ , colSums(is.na(Kona)) < nrow(Kona)]   
                         
Edge <- dist_beta[!is.na(dist_beta$Edge),]
Edge$site <- "Edge"
Edge <- rename(Edge, beta = Edge)
Edge <- Edge[ , colSums(is.na(Edge)) < nrow(Edge)] 
                         
dist_beta <- rbind(Center, Stainbeck, Kona, Edge)
dist_beta <- rename(dist_beta, dist = ï..dist)                         

a <- ggplot(data=dist_beta) + 
  geom_smooth(method='lm', aes(x=logidst, y=beta, colour=site, fill=site, linetype=site), size=1, alpha=0.20)+
  geom_point(aes(x=logidst, y=beta, colour=site), alpha=0.70, size=6, shape=18) + 
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual(values=SiteColors) +
  labs(title="Distance vs zOTU beta diversity", x="Log distance (km)", y="zOTU beta diversity") +
  KipukaTheme +
  guides(color="none", shape="none", fill ="none", linetype="none") +
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
        legend.title = element_text(size=40),
       legend.position = "right")

#################################################################################
#Genetic differentiation wi OTUs vesus distance
dist_diff <- rename(dist_diff, site = ï..Both) 
                         
CC <- dist_diff[!is.na(dist_diff$CenterCenter),]
CC <- rename(CC, diff = CenterCenter)
CC$site<-"Center"                         
CC <- CC[ , colSums(is.na(CC)) < nrow(CC)]                    
                         
EE <- dist_diff[!is.na(dist_diff$EdgeEdge),]
EE <- rename(EE, diff = EdgeEdge)
EE$site<-"Edge"                         
EE <- EE[ , colSums(is.na(EE)) < nrow(EE)]  
                         
HH <- dist_diff[!is.na(dist_diff$HiloHilo),]
HH <- rename(HH, diff = HiloHilo)
HH$site<-"Stainbeck"                         
HH <- HH[ , colSums(is.na(HH)) < nrow(HH)]  
                         
KK <- dist_diff[!is.na(dist_diff$KonaKona),]
KK <- rename(KK, diff = KonaKona)
KK$site<-"Kona"                         
KK <- KK[ , colSums(is.na(KK)) < nrow(KK)] 
                         
LL <- dist_diff[!is.na(dist_diff$Lavalava),]
LL <- rename(LL, diff = Lavalava)
LL$site<-"Lava"                         
LL <- LL[ , colSums(is.na(LL)) < nrow(LL)]                          
                         
dist_diff <- rbind(CC, EE, HH, KK, LL)
                      
b <- ggplot(data=dist_diff) + 
  geom_smooth(method='lm', aes(x=logdist, y=diff, colour=site, fill=site, linetype=site), size=1, alpha=0.20)+
  geom_point(aes(x=logdist, y=diff, colour=site), alpha=0.70, size=6, shape=18) + 
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual(values=SiteColors) +
  labs(title="Distance vs differentiation within OTUs", x="Log distance (km)", y="Differentiation within OTUs") +
  KipukaTheme +
  guides(color="none", shape="none", fill = guide_legend(ncol=1), linetype="none") +
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
        legend.title = element_text(size=40),
       legend.position = "right")

                         
jpeg("Figures/distance.jpg", width=2200, height=1000)   
plot_grid(a, b, ncol = 2, rel_widths = c(1, 1.2))                         
dev.off()
                         
