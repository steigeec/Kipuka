#Project: Kipukas                     #
#Script: NativeNonative_Proportions.R #
#Author: Emma Steigerwald             #
#Date:17 Mar 2022                     #
#######################################

#Set up my working environment
setwd("G:/My Drive/Kipuka")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
font_import()
library(data.table)
library(scales)
library(extrafont)
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

richness <- read.csv("NatNonNat.csv")

#################################################################################
#Proportional representation of nat/nonnat                        
#Stacked bar plot, each site being its own stacked bar, all same height so proportional representation

OTU <- richness[2:nrow(richness),29:814] 
#Make first column the row names
OTU["16",1]<-"ZotuID"
OTU["18",1]<-"3%OTU"
OTU <- OTU[-c(10),]
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)
OTU<-OTU[,-c(2:16)]
i <- c(2:ncol(OTU))                                  
OTU[i] <- apply(OTU[i], 2,            
                    function(x) as.numeric(as.character(x)))
sapply(OTU, class) 

#Now summarize native/non-native species per site using apply
#needs to be data.table for apply to work
OTU<-setDT(OTU)
otu<-OTU[, lapply(.SD, sum), by = INVNAT]
otu<- t(otu)
otu<-as.data.frame(otu)
names(otu) <- c("INV", "NAT")
otu<- otu[-1,]
i <- c(1:ncol(otu))                                  
otu[i] <- apply(otu[i], 2,            
                    function(x) as.numeric(as.character(x)))
otu$p_nat<-otu$NAT/(otu$INV+otu$NAT)
otu$p_non<-otu$INV/(otu$INV+otu$NAT)
otu$Site<-rownames(otu)
otu<-otu[,-c(1:2)]

#Convert wide to longform for plotting
otu<- melt(otu, id.vars=c("Site"))

#Join back other data needed to interpret sites
names(richness)<-richness[c(18),]
richness<-richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))                       

richness<-richness[1:11]
OTU<-merge(richness, otu, by.x="ID", by.y="Site")

#Order as I want panels to appear in plot
OTU$Site <- factor(OTU$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava"))) 
OTU$variable <- factor(OTU$variable, levels = rev(c("p_nat", "p_non")))                 
plotA <-OTU[order(OTU$Site, OTU$Area),]
#Reindex data frame so that it plots this way
rownames(plotA) <- seq(1,nrow(plotA),1)

a<-ggplot(data=plotA, aes(x=reorder(ID,Area), y=value, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + 
  labs(title="A.") +
  facet_grid(cols=vars(Site), scales="free", space="free") +             
  scale_fill_manual("", values=c('#88CCEE', '#44AA99'), labels=c("p_nat"="Native", "p_non"="Invasive"))+
  scale_y_continuous(name="Percent composition", expand = c(0,0), breaks=seq(0,1,.2))+
  #scale_x_discrete(labels=my_labels)+              
  #scale_x_discrete(breaks=OTU$ID, labels=OTU$Area)+ 
  xlab(expression("         [                By increasing size ("~m^2~")                          ]                                                             "))+                             
  KipukaTheme +
  theme(axis.text.y = element_text(angle=45, size=25, vjust=-.5, hjust=.5),
        axis.text.x = element_blank(),
 #       axis.text.x = element_text(size=25, vjust=0.1, hjust=0.1, angle=90),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30), 
        plot.margin = margin(0,0,0,0, "cm"))

                  
A<-ggplot() + 
  geom_boxplot(data=plotA[plotA$variable=="p_non",],aes(x=Site, y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  scale_y_continuous(name="Percent invasive reads")+              
  labs(title="A.", x="") +
  KipukaTheme +
  theme(strip.text = element_text(size = 30), 
        axis.text = element_text(angle=45, size=40), 
        axis.title.y = element_text(size=40), 
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
       legend.position = "right", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))              
                
subsetA<-plotA[plotA$variable=="p_non",]
subsetA<-subsetA[subsetA$Site=="Center" | subsetA$Site== "Edge",]                
A1<-ggplot() + 
  geom_point(data=subsetA,aes(x=as.numeric(Area), y=value, colour=Site), size=6)+
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge")) +
  geom_smooth(method='lm', data=subsetA, aes(x=Area, y=value, colour=Site, fill=Site), size=1, alpha=0.20)+
  scale_fill_manual(values=SiteColors, limits=c("Center", "Edge")) +              
  scale_y_continuous(name="Proportion invasive reads")+              
  labs(title="A.", x="") +
  KipukaTheme +
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +              
  theme(strip.text = element_text(size = 30), 
        axis.text = element_text(angle=45, size=40), 
        axis.title.y = element_text(size=40), 
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
       legend.position = "right", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))    
                     
####################################################################################################
#Also do in terms of species richness of nat/non-nat

richness <- read.csv("NatNonNat.csv")
                
OTU <- richness[2:nrow(richness),29:814] 
#Make first column the row names
OTU["16",1]<-"ZotuID"
OTU["18",1]<-"3%OTU"
OTU <- OTU[-c(10),]
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)
OTU<-OTU[,-c(2:16)]
i <- c(2:ncol(OTU))                                  
OTU[i] <- apply(OTU[i], 2,            
                    function(x) as.numeric(as.character(x)))
sapply(OTU, class) 

#Now summarize native/non-native species per site using apply
#needs to be data.table for apply to work

for (I in 2:ncol(OTU)){
        for (J in 1:nrow(OTU)){
                if (OTU[J,I] > 0) {
                        OTU[J,I]<-1
                }   
                else if(OTU[J,I]==0){
                        OTU[J,I]<-0
               }
               else {
                       print("Something's wrong")
              }               
        }                
}        

OTU<-setDT(OTU)            
otu<-OTU[, lapply(.SD, sum), by = INVNAT]
otu<- t(otu)
otu<-as.data.frame(otu)
names(otu) <- c("INV", "NAT")
otu<- otu[-1,]
i <- c(1:ncol(otu))                                  
otu[i] <- apply(otu[i], 2,            
                    function(x) as.numeric(as.character(x)))
otu$p_nat<-otu$NAT/(otu$INV+otu$NAT)
otu$p_non<-otu$INV/(otu$INV+otu$NAT)
otu$Site<-rownames(otu)
otu<-otu[,-c(1:2)]

#Convert wide to longform for plotting
otu<- melt(otu, id.vars=c("Site"))

#Join back other data needed to interpret sites
names(richness)<-richness[c(18),]
richness<-richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))                       

richness<-richness[1:11]
OTU<-merge(richness, otu, by.x="ID", by.y="Site")

#Order as I want panels to appear in plot
OTU$Site <- factor(OTU$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava"))) 
OTU$variable <- factor(OTU$variable, levels = rev(c("p_nat", "p_non")))                                            
plotB <-OTU[order(OTU$Site, OTU$Area),]
#Reindex data frame so that it plots this way
rownames(plotB) <- seq(1,nrow(plotB),1)                
                
b<-ggplot(data=plotB, aes(x=reorder(ID,Area), y=value, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + 
  labs(title="B.") +
  facet_grid(cols=vars(Site), scales="free", space="free") +             
  scale_fill_manual("", values=c('#88CCEE', '#44AA99'), labels=c("p_nat"="Native", "p_non"="Invasive"))+
  scale_y_continuous(name="Percent composition", expand = c(0,0), breaks=seq(0,1,.2))+
  #scale_x_discrete(labels=my_labels)+              
  #scale_x_discrete(breaks=OTU$ID, labels=OTU$Area)+ 
  xlab(expression("         [                By increasing size ("~m^2~")                          ]                                                             "))+                             
  KipukaTheme +
  theme(axis.text.y = element_text(angle=45, size=25, vjust=-.5, hjust=.5),
        axis.text.x = element_blank(),
 #       axis.text.x = element_text(size=25, vjust=0.1, hjust=0.1, angle=90),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30), 
        plot.margin = margin(0,0,0,0, "cm"))              

B<-ggplot() + 
  geom_boxplot(data=plotB[plotB$variable=="p_non",],aes(x=Site, y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  scale_y_continuous(name="Percent invasive OTUs")+              
  labs(title="B.", x="") +
  guides(fill="none") +             
  KipukaTheme +
  theme(strip.text = element_text(size = 30), 
        axis.text = element_text(angle=45, size=40), 
        axis.title.y = element_text(size=40), 
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

subsetB<-plotB[plotB$variable=="p_non",]
subsetB<-subsetB[subsetB$Site=="Center" | subsetB$Site== "Edge",] 
B1<-ggplot() + 
  geom_point(data=subsetB,aes(x=as.numeric(Area), y=value, colour=Site), size=6)+
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge")) +
  geom_smooth(method='lm', data=subsetB, aes(x=Area, y=value, colour=Site, fill=Site), size=1, alpha=0.20)+
  scale_fill_manual(values=SiteColors, limits=c("Center", "Edge")) +
  scale_y_continuous(name="Proportion invasive OTUs")+              
  labs(title="B.", x="") +
  guides(colour="none", fill="none")+
  KipukaTheme +              
scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +                                    
  theme(strip.text = element_text(size = 30), 
        axis.text = element_text(angle=45, size=40), 
        axis.title.y = element_text(size=40), 
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
        plot.margin = margin(0.2,1,0,1.35, "cm"))
                   
jpeg("Figures/NatNonNat.jpg", width=1500, height=700)
plot_grid(a, b, nrow=2)
dev.off()
                
jpeg("Figures/NatNonNat_box_1.jpg", width=1500, height=1000)
plot_grid(A, B, ncol=2, rel_widths=c(1, .7))
dev.off()                
                
jpeg("Figures/NatNonNat_Scatter.jpg", width=1500, height=1000)
plot_grid(A1, B1, ncol=2, rel_widths=c(1.2, 1))
dev.off()   
 
                   
                   
#################################################################################
#NOW THE SAME FOR ARANEAE ONLY                
#Proportional representation of nat/nonnat                        
#Stacked bar plot, each site being its own stacked bar, all same height so proportional representation
richness <- read.csv("NatNonNat.csv")               
OTU <- richness[2:nrow(richness),29:814] 
#Make first column the row names
OTU["16",1]<-"ZotuID"
OTU["18",1]<-"3%OTU"
OTU <- OTU[-c(10),]
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)
#Filter out Araneae only!!!                
OTU<-OTU[OTU$Order=="Araneae",]               
OTU<-OTU[,-c(2:16)]
i <- c(2:ncol(OTU))                                  
OTU[i] <- apply(OTU[i], 2,            
                    function(x) as.numeric(as.character(x)))
sapply(OTU, class) 

#Now summarize native/non-native species per site using apply
#needs to be data.table for apply to work
OTU<-setDT(OTU)
otu<-OTU[, lapply(.SD, sum), by = INVNAT]
otu<- t(otu)
otu<-as.data.frame(otu)
names(otu) <- c("INV", "NAT")
otu<- otu[-1,]
i <- c(1:ncol(otu))                                  
otu[i] <- apply(otu[i], 2,            
                    function(x) as.numeric(as.character(x)))
otu$p_nat<-otu$NAT/(otu$INV+otu$NAT)
otu$p_non<-otu$INV/(otu$INV+otu$NAT)
otu$Site<-rownames(otu)
otu<-otu[,-c(1:2)]

#Convert wide to longform for plotting
otu<- melt(otu, id.vars=c("Site"))
otu <- otu[otu$Site != "1K08E" & otu$Site != "1K08C",]                      

#Join back other data needed to interpret sites
names(richness)<-richness[c(18),]
richness<-richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))                       

richness<-richness[1:11]
OTU<-merge(richness, otu, by.x="ID", by.y="Site")

#Order as I want panels to appear in plot
OTU$variable <- factor(OTU$variable, levels = rev(c("p_nat", "p_non")))                 
OTU$Area<-round(OTU$Area,10)    
plotA<-OTU                
plotA$Site<-gsub("Stainbeck","Stainback",as.character(plotA$Site)) 
plotA$Site <- factor(plotA$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava")))
               
                
bp <- ggplot() + 
  geom_boxplot(data=plotA[plotA$variable!="p_non",],aes(x=Site, y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  labs(title="A. ", x="", y="Proportion of invasive species") +
  KipukaTheme +
  theme(strip.text = element_text(size = 40), 
        axis.text = element_text(angle=45, size=40), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=40), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=40), 
        legend.title = element_blank(),
       legend.position = "left", 
        plot.margin = margin(0.2,1,0,1.35, "cm"),
       legend.key.height = unit(3, 'cm'), 
        legend.key.width = unit(1, 'cm'))
                
#######################################                
plotA <-OTU[order(OTU$Site, OTU$Area),]
#Reindex data frame so that it plots this way
rownames(plotA) <- seq(1,nrow(plotA),1)
                
a<-ggplot(data=plotA, aes(x=reorder(ID,Area), y=value, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + 
  labs(title="A.") +
  facet_grid(cols=vars(Site), scales="free", space="free") +             
  scale_fill_manual("", values=c('#88CCEE', '#44AA99'), labels=c("p_nat"="Native", "p_non"="Invasive"))+
  scale_y_continuous(name="Percent composition", expand = c(0,0), breaks=seq(0,1,.2))+
  #scale_x_discrete(labels=my_labels)+              
  #scale_x_discrete(breaks=OTU$ID, labels=OTU$Area)+ 
  xlab(expression("         [                By increasing size ("~m^2~")                          ]                                                             "))+                             
  KipukaTheme +
  theme(axis.text.y = element_text(angle=45, size=25, vjust=-.5, hjust=.5),
        axis.text.x = element_blank(),
 #       axis.text.x = element_text(size=25, vjust=0.1, hjust=0.1, angle=90),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30), 
        plot.margin = margin(0,0,0,0, "cm"))

A<-ggplot() + 
  geom_boxplot(data=plotA[plotA$variable=="p_non",],aes(x=reorder(Site, value), y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  scale_y_continuous(name="Percent invasive reads")+              
  labs(title="A.", x="") +
  KipukaTheme +
  theme(strip.text = element_text(size = 30), 
        axis.text = element_text(angle=45, size=40), 
        axis.title.y = element_text(size=40), 
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
       legend.position = "right", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))                    

subsetA<-plotA[plotA$variable=="p_non",]
subsetA<-subsetA[subsetA$Site=="Center" | subsetA$Site== "Edge",]   
subsetA$Site<-gsub("Stainbeck","Stainback",as.character(subsetA$Site))                
A1<-ggplot() + 
  geom_point(data=subsetA,aes(x=as.numeric(Area), y=value, colour=Site), size=6)+
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge")) +
  geom_smooth(method='lm', data=subsetA, aes(x=Area, y=value, colour=Site, fill=Site), size=1, alpha=0.20)+
  scale_fill_manual(values=SiteColors, limits=c("Center", "Edge")) +              
  scale_y_continuous(name="Proportion of invasive reads")+              
  labs(title="B.", x="Kipuka area ("~m^2~")") +
  KipukaTheme +
  guides(colour="none", fill="none")+
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +              
  theme(strip.text = element_text(size = 40), 
        axis.text = element_text(angle=45, size=40), 
        axis.title.y = element_text(size=40), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=40), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=40), 
        legend.title = element_blank(),
       legend.position = "right", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))                  
                
####################################################################################################
#Also do in terms of species richness of nat/non-nat

richness <- read.csv("NatNonNat.csv")
                
OTU <- richness[2:nrow(richness),29:814] 
#Make first column the row names
OTU["16",1]<-"ZotuID"
OTU["18",1]<-"3%OTU"
OTU <- OTU[-c(10),]
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)
#Filter out Araneae only!!!                
OTU<-OTU[OTU$Order=="Araneae",]                  
OTU<-OTU[,-c(2:16)]
i <- c(2:ncol(OTU))                                  
OTU[i] <- apply(OTU[i], 2,            
                    function(x) as.numeric(as.character(x)))
sapply(OTU, class) 

#Now summarize native/non-native species per site using apply
#needs to be data.table for apply to work

for (I in 2:ncol(OTU)){
        for (J in 1:nrow(OTU)){
                if (OTU[J,I] > 0) {
                        OTU[J,I]<-1
                }   
                else if(OTU[J,I]==0){
                        OTU[J,I]<-0
               }
               else {
                       print("Something's wrong")
              }               
        }                
}        

OTU<-setDT(OTU)            
otu<-OTU[, lapply(.SD, sum), by = INVNAT]
otu<- t(otu)
otu<-as.data.frame(otu)
names(otu) <- c("INV", "NAT")
otu<- otu[-1,]
i <- c(1:ncol(otu))                                  
otu[i] <- apply(otu[i], 2,            
                    function(x) as.numeric(as.character(x)))
otu$p_nat<-otu$NAT/(otu$INV+otu$NAT)
otu$p_non<-otu$INV/(otu$INV+otu$NAT)
otu$Site<-rownames(otu)
otu<-otu[,-c(1:2)]

#Convert wide to longform for plotting
otu<- melt(otu, id.vars=c("Site"))
otu <- otu[otu$Site != "1K08E" & otu$Site != "1K08C",]                      

#Join back other data needed to interpret sites
names(richness)<-richness[c(18),]
richness<-richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))                       

richness<-richness[1:11]
OTU<-merge(richness, otu, by.x="ID", by.y="Site")

#Order as I want panels to appear in plot
OTU$Site <- factor(OTU$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava"))) 
OTU$variable <- factor(OTU$variable, levels = rev(c("p_nat", "p_non")))                                 
OTU$Area<-round(OTU$Area,10)                
plotB <-OTU[order(OTU$Site, OTU$Area),]
#Reindex data frame so that it plots this way
rownames(plotB) <- seq(1,nrow(plotB),1)                
                
b<-ggplot(data=plotB, aes(x=reorder(ID,Area), y=value, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + 
  labs(title="B.") +
  facet_grid(cols=vars(Site), scales="free", space="free") +             
  scale_fill_manual("", values=c('#88CCEE', '#44AA99'), labels=c("p_nat"="Native", "p_non"="Invasive"))+
  scale_y_continuous(name="Percent composition", expand = c(0,0), breaks=seq(0,1,.2))+
  #scale_x_discrete(labels=my_labels)+              
  #scale_x_discrete(breaks=OTU$ID, labels=OTU$Area)+ 
  xlab(expression("         [                By increasing size ("~m^2~")                          ]                                                             "))+                             
  KipukaTheme +
  theme(axis.text.y = element_text(angle=45, size=25, vjust=-.5, hjust=.5),
        axis.text.x = element_blank(),
 #       axis.text.x = element_text(size=25, vjust=0.1, hjust=0.1, angle=90),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30), 
        plot.margin = margin(0,0,0,0, "cm"))              
                


B<-ggplot() + 
  geom_boxplot(data=plotB[plotB$variable=="p_non",],aes(x=reorder(Site, value), y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  scale_y_continuous(name="Percent invasive OTUs")+              
  labs(title="B.", x="") +
  guides(fill="none") +             
  KipukaTheme +
  theme(strip.text = element_text(size = 30), 
        axis.text = element_text(angle=45, size=40), 
        axis.title.y = element_text(size=40), 
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

subsetB<-plotB[plotB$variable=="p_non",]
subsetB<-subsetB[subsetB$Site=="Center" | subsetB$Site== "Edge",] 
subsetB$Site<-gsub("Stainbeck","Stainback",as.character(subsetB$Site))                
B1<-ggplot() + 
  geom_point(data=subsetB,aes(x=as.numeric(Area), y=value, colour=Site), size=6)+
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge")) +
  geom_smooth(method='lm', data=subsetB, aes(x=Area, y=value, colour=Site, fill=Site), size=1, alpha=0.20)+
  scale_fill_manual(values=SiteColors, limits=c("Center", "Edge")) +
  scale_y_continuous(name="Proportion of invasive OTUs")+              
  labs(title="C.", x="Kipuka area ("~m^2~")") +
  guides(colour="none", fill="none")+
  KipukaTheme +              
scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +                                    
  theme(strip.text = element_text(size = 40), 
        axis.text = element_text(angle=45, size=40), 
        axis.title.y = element_text(size=40), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=40), 
        plot.title=element_text(size=50),  
        plot.margin = margin(0.2,1,0,1.35, "cm"))
                     
jpeg("Figures/NatNonNat_AraneaeScatter.jpg", width=3000, height=1000)
plot_grid(bp, A1, B1, ncol=3, rel_widths=c(1.15, 1, 1))
dev.off()                   
                
                
jpeg("Figures/NatNonNat_Araneaebox.jpg", width=1500, height=1000)
plot_grid(A, B, ncol=2, rel_widths=c(1, .7))
dev.off()                  
                
jpeg("Figures/NatNonNat_Araneae.jpg", width=1500, height=700)
plot_grid(a, b, nrow=2)
dev.off()             
                
                
#################################################################################
#NOW TRY SPLITTING BY ALL ORDERS............ NOT DONE YET           
#Proportional representation of nat/nonnat                        
#Stacked bar plot, each site being its own stacked bar, all same height so proportional representation
richness <- read.csv("NatNonNat.csv")
                
OTU <- richness[2:nrow(richness),29:814] 
#Make first column the row names
OTU["16",1]<-"ZotuID"
OTU["18",1]<-"3%OTU"
OTU <- OTU[-c(10),]
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)             
OTU<-OTU[,-c(2:16)]
i <- c(2:ncol(OTU))                                  
OTU[i] <- apply(OTU[i], 2,            
                    function(x) as.numeric(as.character(x)))
sapply(OTU, class) 

#Now summarize native/non-native species per site using apply
#needs to be data.table for apply to work
OTU<-setDT(OTU)
otu<-OTU[, lapply(.SD, sum), by = INVNAT]
otu<- t(otu)
otu<-as.data.frame(otu)
names(otu) <- c("INV", "NAT")
otu<- otu[-1,]
i <- c(1:ncol(otu))                                  
otu[i] <- apply(otu[i], 2,            
                    function(x) as.numeric(as.character(x)))
otu$p_nat<-otu$NAT/(otu$INV+otu$NAT)
otu$p_non<-otu$INV/(otu$INV+otu$NAT)
otu$Site<-rownames(otu)
otu<-otu[,-c(1:2)]

#Convert wide to longform for plotting
otu<- melt(otu, id.vars=c("Site"))

#Join back other data needed to interpret sites
names(richness)<-richness[c(18),]
richness<-richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))                       

richness<-richness[1:22]
OTU<-merge(richness, otu, by.x="ID", by.y="Site")

#Order as I want panels to appear in plot
OTU$Site <- factor(OTU$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava"))) 
OTU$variable <- factor(OTU$variable, levels = rev(c("p_nat", "p_non")))                 
OTU$Area<-round(OTU$Area,10)                
plotA <-OTU[order(OTU$Site, OTU$Area),]
#Reindex data frame so that it plots this way
rownames(plotA) <- seq(1,nrow(plotA),1)

a<-ggplot(data=plotA, aes(x=reorder(ID,Area), y=value, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + 
  labs(title="A.") +
  facet_grid(cols=vars(Site), rows=vars(Order), scales="free", space="free") +             
  scale_fill_manual("", values=c('#88CCEE', '#44AA99'), labels=c("p_nat"="Native", "p_non"="Invasive"))+
  scale_y_continuous(name="Percent composition", expand = c(0,0), breaks=seq(0,1,.2))+
  #scale_x_discrete(labels=my_labels)+              
  #scale_x_discrete(breaks=OTU$ID, labels=OTU$Area)+ 
  xlab(expression("         [                By increasing size ("~m^2~")                          ]                                                             "))+                             
  KipukaTheme +
  theme(axis.text.y = element_text(angle=45, size=25, vjust=-.5, hjust=.5),
        axis.text.x = element_blank(),
 #       axis.text.x = element_text(size=25, vjust=0.1, hjust=0.1, angle=90),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30), 
        plot.margin = margin(0,0,0,0, "cm"))
                
####################################################################################################
#Also do in terms of species richness of nat/non-nat

richness <- read.csv("NatNonNat.csv")
                
OTU <- richness[2:nrow(richness),29:814] 
#Make first column the row names
OTU["16",1]<-"ZotuID"
OTU["18",1]<-"3%OTU"
OTU <- OTU[-c(10),]
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)                
OTU<-OTU[,-c(2:16)]
i <- c(2:ncol(OTU))                                  
OTU[i] <- apply(OTU[i], 2,            
                    function(x) as.numeric(as.character(x)))
sapply(OTU, class) 

#Now summarize native/non-native species per site using apply
#needs to be data.table for apply to work

for (I in 2:ncol(OTU)){
        for (J in 1:nrow(OTU)){
                if (OTU[J,I] > 0) {
                        OTU[J,I]<-1
                }   
                else if(OTU[J,I]==0){
                        OTU[J,I]<-0
               }
               else {
                       print("Something's wrong")
              }               
        }                
}        

OTU<-setDT(OTU)            
otu<-OTU[, lapply(.SD, sum), by = INVNAT]
otu<- t(otu)
otu<-as.data.frame(otu)
names(otu) <- c("INV", "NAT")
otu<- otu[-1,]
i <- c(1:ncol(otu))                                  
otu[i] <- apply(otu[i], 2,            
                    function(x) as.numeric(as.character(x)))
otu$p_nat<-otu$NAT/(otu$INV+otu$NAT)
otu$p_non<-otu$INV/(otu$INV+otu$NAT)
otu$Site<-rownames(otu)
otu<-otu[,-c(1:2)]

#Convert wide to longform for plotting
otu<- melt(otu, id.vars=c("Site"))

#Join back other data needed to interpret sites
names(richness)<-richness[c(18),]
richness<-richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))                       

richness<-richness[1:22]
OTU<-merge(richness, otu, by.x="ID", by.y="Site")

#Order as I want panels to appear in plot
OTU$Site <- factor(OTU$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava"))) 
OTU$variable <- factor(OTU$variable, levels = rev(c("p_nat", "p_non")))                                 
OTU$Area<-round(OTU$Area,10)                
plotB <-OTU[order(OTU$Site, OTU$Area),]
#Reindex data frame so that it plots this way
rownames(plotB) <- seq(1,nrow(plotB),1)                
                
b<-ggplot(data=plotB, aes(x=reorder(ID,Area), y=value, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + 
  labs(title="B.") +
  facet_grid(cols=vars(Site), scales="free", space="free") +             
  scale_fill_manual("", values=c('#88CCEE', '#44AA99'), labels=c("p_nat"="Native", "p_non"="Invasive"))+
  scale_y_continuous(name="Percent composition", expand = c(0,0), breaks=seq(0,1,.2))+
  #scale_x_discrete(labels=my_labels)+              
  #scale_x_discrete(breaks=OTU$ID, labels=OTU$Area)+ 
  xlab(expression("         [                By increasing size ("~m^2~")                          ]                                                             "))+                             
  KipukaTheme +
  theme(axis.text.y = element_text(angle=45, size=25, vjust=-.5, hjust=.5),
        axis.text.x = element_blank(),
 #       axis.text.x = element_text(size=25, vjust=0.1, hjust=0.1, angle=90),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30), 
        plot.margin = margin(0,0,0,0, "cm"))              
                
jpeg("Figures/NatNonNat_Araneae.jpg", width=1500, height=1000)
plot_grid(a, b, nrow=2)
dev.off()                
