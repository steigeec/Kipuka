#Project: Kipukas                 #
#Script: Proportional_Represent.R #
#Author: Emma Steigerwald         #
#Date:17 Feb 2022                 #
###################################

#Set up my working environment
setwd("G:/My Drive/Kipuka")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
font_import()


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
#Proportional representation of orders by site                         
#Stacked bar plot, each site being its own stacked bar, all same height so proportional representation

rep <- richness[1:22]
i<-seq(15, 21, 1)
rep[ , i] <- apply(rep[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))                         
rep$totalRichness <- rep$Hemiptera + rep$Lepidoptera + rep$Pscoptera +  rep$Araneae + rep$Coleoptera + rep$Diptera
                   
rep <- melt(rep, idvars = c("SiteID", "Site", "totalRichness"), measure.vars = c("Araneae", "Pscoptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Acari", "Coleoptera", "Diptera"))
rep$my_site <- paste(rep$Site, rep$SiteID)
rep$prop <- rep$value/rep$totalRichness
                   
#I want ordered by my sites
rep$Site <- factor(rep$Site, levels = rev(c("Kona","Stainbeck",  "Center", "Edge", "Lava")))  
rep <- rename(rep, id = ï..ID) 
rep$Arealog<-round(as.numeric(rep$Arealog),1)
rep<-rep %>% dplyr::mutate(Arealog = tidyr::replace_na(Arealog, ""))                   
                
a<- ggplot(data=rep, aes(x=reorder(my_site, Area), y=prop, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + #, size=0.4, key_glyph = "polygon3"
  labs(title="A.") +
  xlab(expression("         [                 By increasing size ("~m^2~")                         ]                                                             "))+                 
  facet_grid(cols=vars(Site), rows=vars(variable), scales="free", space="free") +             
  scale_fill_manual("Taxon", values=c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', "black"))+
  scale_y_continuous(name="proportional zOTU representation", expand = c(0,0.1), breaks=seq(0,.60,.20))+
  #geom_text(aes(label = round(Arealog,1)),vjust=-.25, size=8) +
  #scale_x_discrete(breaks=rep$my_site, labels=rep$Arealog)+                 
  KipukaTheme +
  theme(axis.text.x = element_blank(), #element_text(angle=45, size=25, vjust=-.1, hjust=1),
        axis.text.y = element_text(size=25, vjust=0.5, hjust=0.5),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30))

b<- ggplot(data=rep, aes(x=reorder(my_site, Area), y=value, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + #, size=0.4, key_glyph = "polygon3"
  labs(title="B.") +
  xlab(expression("         [                 By increasing size ("~m^2~")                         ]                                                             "))+                 
  facet_grid(cols=vars(Site), rows=vars(variable), scales="free", space="free") +             
  scale_fill_manual("Taxon", values=c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', "black"))+
  scale_y_continuous(name="zOTU count per order", expand = c(0,0.1), breaks=seq(0,100,20))+
  #geom_text(aes(label = round(Arealog,1)),vjust=-.25, size=8) +
  #scale_x_discrete(breaks=rep$my_site, labels=rep$Arealog)+                 
  KipukaTheme +
  theme(axis.text.x = element_blank(), #element_text(angle=45, size=25, vjust=-.1, hjust=1),
        axis.text.y = element_text(size=25, vjust=0.5, hjust=0.5),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30))

                   
jpeg("Figures/Order_Representation.jpg", width=1500, height=2000)
plot_grid(a, b, nrow=2, rel_heights=c(1, .75))                   
dev.off()
