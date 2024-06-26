
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

dist_beta_0 <- read.csv("Distance_v_beta.csv")
otu <- read.csv("Distances_Without_Kona8andsmall_kipuka.csv")
dist_diff <- read.csv("Distance_v_differentiation.csv")
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
colnames(acari_beta) <- c("Geo", "log.1", "Kona", "Stainbeck", "Center", "Edge")
acari_beta<-as.data.frame(acari_beta[,c(1, 3, 4, 5, 6)])

#Now rearrange from wide to long format
library(reshape2)
acari_beta <- melt(acari_beta, id= c("Geo"))

#Assign to "kipuka" or "forest"
acari_beta$group[acari_beta$variable=="Center" | acari_beta$variable=="Edge"] <- "Kipuka"
acari_beta$group[acari_beta$variable=="Kona" | acari_beta$variable=="Stainbeck"] <- "Continuous forest"

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

b <- ggplot(data=dist_beta) + 
  geom_smooth(method='lm', aes(x=logidst, y=beta, colour=site, fill=site, linetype=site), size=1, alpha=0.20)+
  geom_point(aes(x=logidst, y=beta, colour=site), alpha=0.70, size=6, shape=18) + 
  facet_wrap(~site, ncol=4)+
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
################################################################################                
                         
jpeg("Figures/Fig_4.jpg", width=2000, height=1000)   
a
#plot_grid(a, b, nrow = 2, rel_heights = c(1, 1))                         
dev.off()
                                    
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
                         
dist_diff <- rbind(CC, EE, HH, KK) 
                      
c <- ggplot(data=dist_diff) + 
  geom_smooth(method='lm', aes(x=logdist, y=diff, colour=site, fill=site), size=1, alpha=0.20)+ #, linetype=site
  geom_point(aes(x=logdist, y=diff, colour=site), alpha=0.70, size=6, shape=18) + 
  facet_wrap(~site, ncol=4) +                       
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual(values=SiteColors) +
  labs(title="C.", x="Distance (km)", y="Haplotype differentiation") +
  KipukaTheme +
  guides(color="none", shape="none", fill="none") +
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

#Now summarize beta diversity between site types
acari_beta<-acari_beta[,c(4, 5, 7)]
acari_beta$metric <- "3% OTU"
dist_beta<-dist_beta[,c(2, 3, 4)]
dist_beta$metric <- "zOTU"
names(dist_beta)<-c( "log_dist", "dist", "Site.x", "metric")
dist_diff<-dist_diff[,c(1, 3, 4)] 
dist_diff$metric <- "Haplotype"
names(dist_diff)<-c("Site.x", "log_dist", "dist",  "metric")

beta <- rbind(dist_diff, dist_beta, acari_beta) 
#Reorder facets
beta$metric <- factor(beta$metric, levels = rev(c("Haplotype", "zOTU", "3% OTU")))                     
                   
jpeg("Figures/BetaSummary.jpg", width=1000, height=1000)                   
ggplot() + 
  geom_boxplot(data=beta,aes(x=reorder(Site.x, dist), y=dist, fill=Site.x), color="black", size=1)+
  facet_wrap(~metric, scales="free") +
  scale_fill_manual(values=SiteColors) +
  labs(title="Beta diversity by site type", x="") +
  KipukaTheme +
  theme(strip.text = element_text(size = 30), 
        axis.text = element_text(angle=45, size=40), 
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
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))
dev.off()
