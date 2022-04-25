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
library(scales)
font_import()


#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainback"="#999933", "C-E"="black", "C-F"="white")
#Establish some themes up top to apply to all
KipukaTheme <- theme(axis.title=element_text(size=50), 
        axis.text = element_text(size=25, angle=50), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=40), 
        legend.key.height = unit(1, "cm"), 
        legend.key.width = unit(1.5,"cm"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.title = element_text(size=40), 
        text = element_text(family = "serif"), 
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.spacing.y = unit(0.1,"cm")) 

                   
#################################################################################
#Version with all the data

richness <- read.csv("merged_by_site_2.csv")
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))
dist_beta_0 <- read.csv("Distance_v_beta.csv")
otu <- read.csv("3OTU.csv")
dist_diff <- read.csv("Distance_v_differentiation.csv")
geo_dist<-read.csv("geo_dist.csv")
nmds<-read.csv("nmds3otu.csv")
#Fix spelling error on sheet before proceeding
nmds$Site<-gsub("Stainbeck","Stainback",as.character(nmds$Site))
BC<-read.csv("BC.csv")
BC3<-read.csv("BC3.csv")
CE<-read.csv("C-E.csv")

#NMDS plot for 3%OTU

#Create an NMDS plot with columns MDS1 and MDS2
nmds$Area<-as.numeric(gsub(",","",as.character(nmds$Area)))

nmds$pointsize<-round(sqrt(nmds$Area)/10,0)
nmds$pointsize[nmds$Site=="Lava" & is.na(nmds$pointsize)] <- 2
nmds$pointsize[nmds$Site=="Kona" & is.na(nmds$pointsize)] <- 2
nmds$pointsize[nmds$Site=="Stainback" & is.na(nmds$pointsize)] <- 2

nmds <- nmds[nmds$Site != "1K08E" & nmds$Site != "1K08C",] 

a <- ggplot() + 
  geom_point(data=nmds[nmds$Site!=c("C-E", "C-F"),],aes(x=MDS1OTU,y=MDS2OTU,colour=Site, size=pointsize, shape=Site), alpha=0.70, stroke=3) + 
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge", "Lava", "Kona", "Stainback")) +
  scale_shape_manual("Site", values=c("Center" = 16, "Edge" = 16, "Lava"=3, "Kona"=2, "Stainback"=2)) +
  scale_size_continuous("Kipuka area ("~m^2~")", range=c(2,32), breaks=seq(2,32,5), labels=round((10*seq(2,32,5))^2,100)) +
  labs(title="A.", x="NMDS1", y="NMDS2") +
  #coord_equal() +
  scale_x_continuous(breaks=seq(-2,1.5,0.5)) +
  guides(colour = guide_legend(override.aes = list(size=4))) + 
  KipukaTheme +
  theme(panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5),
        plot.title=element_text(size=50))

######################################################33
#Beta diversity summary plot (B)

#Grab beta diversity between centers and each continuous forest type
#First, zOTU between them
colnames(BC)<- gsub('[X]', '', colnames(BC))
rownames(BC)<-BC[,1]
BC<-BC[,-1]
df <- data.frame(matrix(ncol = 4, nrow = 0))
names(df)<-c("log_dist", "dist", "Site.x", "metric")
out_row<-1
for (ROW in 1:nrow(BC)){
        for (COL in 1:ncol(BC)){             
                df[out_row,2]<-BC[ROW,COL] 
                df[out_row,4]<-"zOTU"
                df[out_row,3]<-"C-F"
                out_row<-out_row+1
        }        
}

#Now for 3%OTU
colnames(BC3)<- gsub('[X]', '', colnames(BC3))
rownames(BC3)<-BC3[,1]
BC3<-BC3[,-1]
df3 <- data.frame(matrix(ncol = 4, nrow = 0))
names(df3)<-c("log_dist", "dist", "Site.x", "metric")
out_row<-1
for (ROW in 1:nrow(BC3)){
        for (COL in 1:ncol(BC3)){             
                df3[out_row,2]<-BC3[ROW,COL] 
                df3[out_row,4]<-"3% OTU"
                df3[out_row,3]<-"C-F"
                out_row<-out_row+1
        }        
}

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
#Remove all the XXXX from colnames
colnames(acari_beta)<- gsub('[X]', '', colnames(acari_beta))
rownames(acari_beta)<-acari_beta[,1]
acari_beta<-acari_beta[,-1]
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
#Remove flipped pairs
acari_beta <- acari_beta%>% distinct(dist, .keep_all= TRUE)
acari_beta$dist<-as.numeric(acari_beta$dist)

#Maybe join together tables of each part...?
dist_beta <- dist_beta_0[,1:7]
i<-c(3, 4, 5, 6, 7)
dist_beta[ , i] <- apply(dist_beta[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))

Center <- dist_beta[!is.na(dist_beta$Center),]
Center$site <- "Center"
Center <- rename(Center, c("beta"="Center"))
Center <- Center[ , colSums(is.na(Center)) < nrow(Center)]                    
                         
Stainbeck <- dist_beta[!is.na(dist_beta$Stainbeck),]
Stainbeck$site <- "Stainbeck"  
Stainbeck <- rename(Stainbeck, c("beta"="Stainbeck"))
Stainbeck <- Stainbeck[ , colSums(is.na(Stainbeck)) < nrow(Stainbeck)]   
                         
Kona <- dist_beta[!is.na(dist_beta$Kona),]
Kona$site <- "Kona"
Kona <- rename(Kona, c("beta"="Kona"))
Kona <- Kona[ , colSums(is.na(Kona)) < nrow(Kona)]   
                         
Edge <- dist_beta[!is.na(dist_beta$Edge),]
Edge$site <- "Edge"
Edge <- rename(Edge, c("beta"="Edge"))
Edge <- Edge[ , colSums(is.na(Edge)) < nrow(Edge)] 
                         
Lava <- dist_beta[!is.na(dist_beta$Lava),]
Lava$site <- "Lava"
Lava <- rename(Lava, c("beta"="Lava"))
Lava <- Lava[ , colSums(is.na(Lava)) < nrow(Lava)]                          
                         
dist_beta <- rbind(Center, Stainbeck, Kona, Edge, Lava)
dist_beta <- rename(dist_beta, c("dist"="ï..dist"))                         


#Now summarize beta diversity between site types
acari_beta<-acari_beta[,c(4, 5, 7)]
acari_beta$metric <- "3% OTU"
dist_beta<-dist_beta[,c(2, 3, 4)]
dist_beta$metric <- "zOTU"
names(dist_beta)<-c("dist", "log_dist", "Site.x", "metric")

beta <- rbind(dist_beta, acari_beta) 
#Reorder facets
beta$metric <- factor(beta$metric, levels = rev(c("zOTU", "3% OTU")))   
beta <- beta[beta$Site != "1K08E" & beta$Site != "1K08C",]
beta <- beta[beta$Site.x != "C-F" & beta$Site.x != "C-E",]

#Fix spelling error on sheet before proceeding
beta$Site.x<-gsub("Stainbeck","Stainback",as.character(beta$Site.x))
                         
beta$Site.x<-as.factor(beta$Site.x)
beta$Site.x <- factor(beta$Site.x, levels=c("Lava", "Edge", "Center", "Stainback", "Kona"))  
                         
                         
b<- ggplot() + 
  geom_boxplot(data=beta[beta$metric=="3% OTU" & !is.na(beta$Site.x),],aes(x=Site.x, y=dist, fill=Site.x), color="black", size=1)+
  #facet_wrap(~metric, scales="free") +
  scale_fill_manual(values=SiteColors) +
  labs(title="B.", y="3% OTU beta diversity", x="") +
  KipukaTheme +
  guides(fill="none")+                       
  theme(strip.text = element_text(size = 50), 
        axis.text.x = element_text(angle=45, size=50, hjust=1, vjust=1), 
        axis.text.y = element_text(angle=45, size=50), 
        axis.title.y = element_text(size = 50),
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
        legend.text=element_text(size=40), 
        legend.title = element_blank(),
       legend.position = "top", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))

#######################################################
#Bray curtis center-edge pairs only
CE<-merge(CE, richness, by.x="log_dist", by.y="ï..ID")

#Fix spelling error on sheet before proceeding
CE$Site<-gsub("Stainbeck","Stainback",as.character(CE$Site))                         
                         
c<-ggplot(data=CE[CE$metric=="3% OTU",]) + 
  geom_smooth(method='lm', aes(x=Area, y=dist), colour="black", size=1, alpha=0.20)+
  geom_point(aes(x=Area, y=dist), colour="#6699CC", fill="#332288", alpha=0.70, size=8, shape=21, stroke=7) + 
  labs(title="C.", x="Kipuka area ("~m^2~")", y="3% OTU beta diversity") +
  KipukaTheme +
  #coord_cartesian(ylim=c(0.4, 1))+
  guides(colour="none")+
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +                   
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


#######################################################
jpeg("Figures/NMDS-turnover.jpg", width=3000, height=1000) 
plot_grid(a,b,c, nrow=1, rel_widths=c(1.3, .9, 1), rel_heights=c(1,1,1))                         
dev.off()

                     
                     
                     
                     
                     
#################################################################################
#Version without the 8 smallest Kipuka

richness <- read.csv("merged_by_site_2.csv")
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))
dist_beta_0 <- read.csv("Distance_v_beta.csv")
otu <- read.csv("3OTU.csv")
dist_diff <- read.csv("Distance_v_differentiation.csv")
geo_dist<-read.csv("geo_dist.csv")
nmds<-read.csv("nmds3otu.csv")
BC<-read.csv("BC.csv")
BC3<-read.csv("BC3.csv")
CE<-read.csv("C-E.csv")

#NMDS plot for 3%OTU

#Create an NMDS plot with columns MDS1 and MDS2
nmds$Area<-as.numeric(gsub(",","",as.character(nmds$Area)))
nmds$pointsize<-round(sqrt(nmds$Area)/10,0)
nmds$pointsize[nmds$Site=="Lava" & is.na(nmds$pointsize)] <- 2
nmds$pointsize[nmds$Site=="Kona" & is.na(nmds$pointsize)] <- 2
nmds$pointsize[nmds$Site=="Stainbeck" & is.na(nmds$pointsize)] <- 2

nmds <- nmds[nmds$ID != "1K08E" & nmds$ID != "1K08C",] 
nmds <- nmds[nmds$Site=="Lava" | nmds$Site=="Kona" | nmds$Site=="Stainbeck" |  nmds$Area>5000,]                      

a <- ggplot() + 
  geom_point(data=nmds[nmds$Site!=c("C-E", "C-F"),],aes(x=MDS1OTU,y=MDS2OTU,colour=Site, size=pointsize, shape=Site), alpha=0.70, stroke=3) + 
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge", "Lava", "Kona", "Stainbeck")) +
  scale_shape_manual("Site", values=c("Center" = 16, "Edge" = 16, "Lava"=3, "Kona"=2, "Stainbeck"=2)) +
  scale_size_continuous("Kipuka area ("~m^2~")", range=c(2,32), breaks=seq(2,32,5), labels=round((10*seq(2,32,5))^2,100)) +
  labs(title="A.", x="NMDS1", y="NMDS2") +
  #coord_equal() +
  scale_x_continuous(breaks=seq(-2,1.5,0.5)) +
  guides(colour = guide_legend(override.aes = list(size=4))) + 
  KipukaTheme +
  theme(panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5),
        plot.title=element_text(size=50))

######################################################
#Beta diversity summary plot (B)

#Grab beta diversity between centers and each continuous forest type
#First, zOTU between them
colnames(BC)<- gsub('[X]', '', colnames(BC))
rownames(BC)<-BC[,1]
BC<-BC[,-1]
df <- data.frame(matrix(ncol = 4, nrow = 0))
names(df)<-c("log_dist", "dist", "Site.x", "metric")
out_row<-1
for (ROW in 1:nrow(BC)){
        for (COL in 1:ncol(BC)){             
                df[out_row,2]<-BC[ROW,COL] 
                df[out_row,4]<-"zOTU"
                df[out_row,3]<-"C-F"
                out_row<-out_row+1
        }        
}

#Now for 3%OTU
colnames(BC3)<- gsub('[X]', '', colnames(BC3))
rownames(BC3)<-BC3[,1]
BC3<-BC3[,-1]
df3 <- data.frame(matrix(ncol = 4, nrow = 0))
names(df3)<-c("log_dist", "dist", "Site.x", "metric")
out_row<-1
for (ROW in 1:nrow(BC3)){
        for (COL in 1:ncol(BC3)){             
                df3[out_row,2]<-BC3[ROW,COL] 
                df3[out_row,4]<-"3% OTU"
                df3[out_row,3]<-"C-F"
                out_row<-out_row+1
        }        
}

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
#Remove all the XXXX from colnames
colnames(acari_beta)<- gsub('[X]', '', colnames(acari_beta))
rownames(acari_beta)<-acari_beta[,1]
acari_beta<-acari_beta[,-1]
acari_beta<-data.frame(col=colnames(acari_beta)[col(acari_beta)], row=rownames(acari_beta)[row(acari_beta)], dist=c(acari_beta))
#Add attributes of each site
#First we add whether it's center, edge, etc etc etc
acari_beta<-merge(acari_beta, richness[,c(1,9,10)], by.x="col", by.y="ï..ID")
acari_beta<-merge(acari_beta, richness[,c(1,9)], by.x="row", by.y="ï..ID")
#If Site.x and Site.y are not the same (e.g. not both Center and Center), then throw out that row
acari_beta<-acari_beta[acari_beta$Site.x==acari_beta$Site.y,]

#Now add the distances between these sites
acari_beta$index<-paste(acari_beta$row, acari_beta$col, sep="_")
acari_beta<-merge(acari_beta, geo_dist, by.x="index", by.y="index")
#remove the same-site pairs
acari_beta<-acari_beta[acari_beta$row!=acari_beta$col,]     
#Remove flipped pairs
acari_beta <- acari_beta%>% distinct(dist, .keep_all= TRUE)
acari_beta$dist<-as.numeric(acari_beta$dist)
acari_beta<-acari_beta[acari_beta$Site.x=="Lava" | acari_beta$Site.x=="Kona" | acari_beta$Site.x=="Stainbeck" |  acari_beta$Area>5000,]                     

#Maybe join together tables of each part...?
dist_beta <- dist_beta_0[,1:7]
i<-c(3, 4, 5, 6, 7)
dist_beta[ , i] <- apply(dist_beta[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))

Center <- dist_beta[!is.na(dist_beta$Center),]
Center$site <- "Center"
Center <- rename(Center, c("beta"="Center"))
Center <- Center[ , colSums(is.na(Center)) < nrow(Center)]                    
                         
Stainbeck <- dist_beta[!is.na(dist_beta$Stainbeck),]
Stainbeck$site <- "Stainbeck"  
Stainbeck <- rename(Stainbeck, c("beta"="Stainbeck"))
Stainbeck <- Stainbeck[ , colSums(is.na(Stainbeck)) < nrow(Stainbeck)]   
                         
Kona <- dist_beta[!is.na(dist_beta$Kona),]
Kona$site <- "Kona"
Kona <- rename(Kona, c("beta"="Kona"))
Kona <- Kona[ , colSums(is.na(Kona)) < nrow(Kona)]   
                         
Edge <- dist_beta[!is.na(dist_beta$Edge),]
Edge$site <- "Edge"
Edge <- rename(Edge, c("beta"="Edge"))
Edge <- Edge[ , colSums(is.na(Edge)) < nrow(Edge)] 
                         
Lava <- dist_beta[!is.na(dist_beta$Lava),]
Lava$site <- "Lava"
Lava <- rename(Lava, c("beta"="Lava"))
Lava <- Lava[ , colSums(is.na(Lava)) < nrow(Lava)]                          
                         
dist_beta <- rbind(Center, Stainbeck, Kona, Edge, Lava)
dist_beta <- rename(dist_beta, c("dist"="ï..dist"))                         


#Now summarize beta diversity between site types
acari_beta<-acari_beta[,c(4, 5, 8)]
acari_beta$metric <- "3% OTU"
dist_beta<-dist_beta[,c(2, 3, 4)]
dist_beta$metric <- "zOTU"
names(dist_beta)<-c( "log_dist", "dist","Site.x", "metric")

beta <- rbind(dist_beta, acari_beta, CE, df, df3) 
#Reorder facets
beta$metric <- factor(beta$metric, levels = rev(c("zOTU", "3% OTU")))   
beta <- beta[beta$Site != "1K08E" & beta$Site != "1K08C",]
beta <- beta[beta$Site.x != "C-F" & beta$Site.x != "C-E",]

beta$Site.x<-as.factor(beta$Site.x)
beta$Site.x <- factor(beta$Site.x, levels=c("Lava", "Edge", "Center", "Stainbeck", "Kona"))                         
                         
b<- ggplot() + 
  geom_boxplot(data=beta[beta$metric=="3% OTU",],aes(x=Site.x, y=dist, fill=Site.x), color="black", size=1)+
  #facet_wrap(~metric, scales="free") +
  scale_fill_manual(values=SiteColors) +
  labs(title="B.", x="") +
  KipukaTheme +
  guides(fill="none")+                       
  theme(strip.text = element_text(size = 50), 
        axis.text.x = element_text(angle=45, size=50, hjust=1, vjust=1), 
        axis.text.y = element_text(angle=45, size=50), 
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
        legend.text=element_text(size=40), 
        legend.title = element_blank(),
       legend.position = "top", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))

#######################################################
#Bray curtis center-edge pairs only
CE<-merge(CE, richness, by.x="log_dist", by.y="ï..ID")
CE<- CE[CE$Area>5000,]                                              
                         
c<-ggplot(data=CE[CE$metric=="3% OTU",]) + 
  geom_smooth(method='lm', aes(x=Area, y=dist), colour="black", size=1, alpha=0.20)+
  geom_point(aes(x=Area, y=dist), colour="#6699CC", fill="#332288", alpha=0.70, size=8, shape=21, stroke=5) + 
  labs(title="C.", x="Kipuka area ("~m^2~")", y="3% OTU beta diversity") +
  KipukaTheme +
  #coord_cartesian(ylim=c(0.4, 1))+
  guides(colour="none")+
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +                   
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


#######################################################
jpeg("Figures/NMDS-turnover-excludingsmallest8.jpg", width=3000, height=1000) 
plot_grid(a,b,c, nrow=1, rel_widths=c(1.3, .9, 1), rel_heights=c(1,1,1))                         
dev.off()













#################################################################################
#NMDS plot for zOTU

#Create an NMDS plot with columns MDS1 and MDS2
richness_mod <- richness
richness_mod$Area<-as.numeric(gsub(",","",as.character(richness_mod$Area)))
richness_mod$pointsize<-round(sqrt(richness_mod$Area)/10,0)
richness_mod$pointsize[richness_mod$Site=="Lava" & is.na(richness_mod$pointsize)] <- 2
richness_mod$pointsize[richness_mod$Site=="Kona" & is.na(richness_mod$pointsize)] <- 2
richness_mod$pointsize[richness_mod$Site=="Stainbeck" & is.na(richness_mod$pointsize)] <- 2
richness_mod <- richness_mod[richness_mod$Site != "1K08E" & richness_mod$Site != "1K08C",] 

jpeg("Figures/NMDS_zOTU.jpg", width=1400, height=1100)
ggplot() + 
  geom_point(data=richness_mod,aes(x=MDS1,y=MDS2,colour=Site, size=pointsize, shape=Site), alpha=0.70, stroke=3) + 
  scale_colour_manual(values=SiteColors) +
  scale_shape_manual("Site", values=c("Center" = 16, "Edge" = 16, "Lava"=3, "Kona"=2, "Stainbeck"=2)) +
  scale_size_continuous("Kipuka area ("~m^2~")", range=c(2,32), breaks=seq(2,32,5), labels=round((10*seq(2,32,5))^2,100)) +
  labs(title="", x="NMDS1", y="NMDS2") +
  #coord_equal() +
  scale_y_continuous(limits=c(-1,1.20)) +
  #guides(colour = "none", size="none", shape="none") + 
  KipukaTheme +
  theme(panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5),
        plot.title=element_text(size=50))
dev.off()

