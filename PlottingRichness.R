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

#Paste all these various dataframes together
order_all<-rbind(acari_beta, Araneae_beta, Coleoptera_beta, Diptera_beta, Hemiptera_beta, Lepidoptera_beta, Psocoptera_beta) 


jpeg("Figures/Order_beta_diversity.jpg", width=1500, height=2000)
ggplot(data=order_all) + 
  geom_smooth(method='lm', aes(x=log_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20)+
  geom_point(aes(x=log_dist, y=dist, colour=Site.x), alpha=0.70, size=4, shape=18) + 
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual(values=SiteColors) +
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
  scale_fill_manual(values=SiteColors) +
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
  labs(title="Distance vs zOTU beta diversity", x="Distance (km)", y="zOTU beta diversity") +
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
  labs(title="Distance vs differentiation within OTUs", x="Distance (km)", y="Differentiation within OTUs") +
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

jpeg("Figures/Order_Richness.jpg", width=4000, height=2000)
ggplot() + 
  geom_boxplot(data=richness_mod_0,aes(x=Arealog, y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  facet_wrap(~variable, nrow=2, scales="free") +
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
#Proportional representation of orders by site                         
#Stacked bar plot, each site being its own stacked bar, all same height so proportional representation

rep <- richness[1:22]
i<-seq(15, 21, 1)
rep[ , i] <- apply(rep[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))                         
rep$totalRichness <- rep$Hemiptera + rep$Hymenoptera + rep$Lepidoptera + rep$Pscoptera +  rep$Acari +  rep$Araneae + rep$Coleoptera + rep$Diptera
                   
rep <- melt(rep, idvars = c("SiteID", "Site", "totalRichness"), measure.vars = c("Araneae", "Pscoptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Acari", "Coleoptera", "Diptera"))
rep$my_site <- paste(rep$Site, rep$SiteID)
rep$prop <- rep$value/rep$totalRichness
                   
#I want ordered by my sites
rep$Site <- factor(rep$Site, levels = rev(c("Kona","Stainbeck",  "Center", "Edge", "Lava")))  
rep <- rename(rep, id = ï..ID)                   
                
jpeg("Figures/Order_Representation.jpg", width=1500, height=1000)
ggplot(data=rep, aes(x=reorder(my_site, Arealog), y=prop, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + #, size=0.4, key_glyph = "polygon3"
  labs(title="Proportion") +
  facet_grid(cols=vars(Site), scales="free", space="free") +             
  scale_fill_manual(values=c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', "black"))+
  scale_y_continuous(name="Proportion", limits=c(0,1.01), expand = c(0,0))+
  KipukaTheme +
  theme(axis.text = element_text(angle=45, size=15, vjust=-1, hjust=1), strip.text = element_text(size = 30))
dev.off()

#################################################################################
#Kipuka size versus SR/SROTU

richness_mod <- richness[richness$Site=="Center" | richness$Site=="Edge",]
richness_mod <- melt(richness_mod, idvars = c("SiteID", "Arealog"), measure = c("SR", "SROTU"))


a <- ggplot() + 
  geom_smooth(method='lm', data=richness_mod,aes(x=Arealog, y=value, colour=Site, linetype=variable, fill=Site), size=1, alpha=0.20)+
  geom_point(data=richness_mod,aes(x=Arealog, y=value, colour=Site, shape=variable), alpha=0.70, size=6, stroke = 3) + 
  scale_shape_manual("Site", values=c("SR" = 0, "SROTU"=15)) +
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual(values=SiteColors)+                 
  #scale_linetype_discrete(values=c(2,5)) +
  labs(title="Size vs species richness", x="Log area ("~km^2~")", y="Species richness") +
  KipukaTheme +
  guides(color="none", shape="none", fill="none") +
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
        legend.title = element_blank(),
       legend.position = "top")

################################################################################
#haplotype richness within OTUs for Kipuka Centers and Kipuka edges

richness_mod_1 <- richness[richness$Site=="Center" | richness$Site=="Edge",]

b <- ggplot() + 
  geom_smooth(method='lm', data=richness_mod_1,aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site, fill=Site), size=1, alpha=0.20)+
  geom_point(data=richness_mod_1,aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site), shape=18, alpha=0.70, size=6, stroke = 3) + 
  facet_wrap(~Site)+
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual(values=SiteColors, limits = c("Center", "Edge"))+                 
  labs(title="Size vs haplotype richness within OTUs", x="Log area ("~km^2~")", y="Haplotype richness within OTUs") +
  KipukaTheme +
  guides(color="none") +
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
        legend.title = element_blank(),
       legend.position = "top")

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
        axis.text = element_text(angle=45, size=40), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=45), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_blank(),
       legend.position = "top")

###########################################################
#Plot these three together

jpeg("Figures/Figure2.jpg", width=3000, height=1000)
plot_grid(a, b, c, ncol = 3, rel_widths = c(1, 1, 1.5))
dev.off()

#################################################################################
#NMDS plot

#Create an NMDS plot with columns MDS1 and MDS2
richness_mod <- richness
richness_mod$Arealog[richness_mod$Site=="Lava" & is.na(richness_mod$Arealog)] <- 2
richness_mod$Arealog[richness_mod$Site=="Kona" & is.na(richness_mod$Arealog)] <- 9
richness_mod$Arealog[richness_mod$Site=="Stainbeck" & is.na(richness_mod$Arealog)] <- 9
richness_mod$Arealog<-richness_mod$Arealog^1.5

jpeg("Figures/NMDS_1.jpg", width=1000, height=1000)
ggplot() + 
  geom_point(data=richness_mod,aes(x=MDS1,y=MDS2,colour=Site, size=(Arealog), shape=Site), alpha=0.70, stroke=3) + 
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

