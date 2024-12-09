#Set up my working environment
setwd("H:/My Drive/Kipuka/Data")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
library(scales)
library(mgcv)
#install.packages("remotes")
#remotes::install_github("Jtrachsel/funfuns")
library(funfuns)
library(BiocManager)
library(vegan)
#BiocManager::install("phyloseq")
library(phyloseq)
library(tidyr)

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainback"="#999933", "C-E"="black", "C-F"="white")
  
#Establish some themes up top to apply to all
KipukaTheme <- theme(axis.title=element_text(size=70), 
        axis.text = element_text(size=70, angle=50), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        plot.title=element_text(size=70), 
        legend.text=element_text(size=70), 
        legend.key.height = unit(1, "cm"), 
        legend.key.width = unit(1.5,"cm"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.title = element_text(size=70), 
        text = element_text(family = "serif"), 
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.spacing.y = unit(0.1,"cm")) 

                   
#################################################################################
# Read in required data... 
richness <- read.csv("merged_by_site_2.csv")
names(richness) <- richness[18,]
richness <- richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))
richness$Site<-gsub("Stainbeck","Stainback",as.character(richness$Site))
# beta diversity for 3 % OTU for within-area
otu <- read.csv("OTU3_Bray_Order.csv")
# beta diversity for zOTU for within-area
zOTUbeta <- read.csv("zOTU_jaccard_Order.csv")
# geographic distances
geo_dist<-read.csv("geo_dist.csv")

#######

#Wrangel geographic distances
rownames(geo_dist) <- geo_dist[,1]
geo_dist <- geo_dist[,-1]
#remove x from all the column names
names(geo_dist)<-sub("X*", "", names(geo_dist))
geo_dist<-as.matrix(geo_dist)
# Turn matrix to long format
geo_dist<-data.frame(col=colnames(geo_dist)[col(geo_dist)], row=rownames(geo_dist)[row(geo_dist)], dist=c(geo_dist))
#make an index column that reps this particular combination of sites
geo_dist$index<-paste(geo_dist$col, geo_dist$row, sep="_")
#provide a log-transformed geo distance
geo_dist$log_dist<-log(geo_dist$dist+0.00001)
# geo_dist<-geo_dist[,c(4, 5)]

###########

# merge zOTU and 3OTU data
otu$index <- paste(otu$Var1, otu$Var2, sep="_")
zOTUbeta$index <- paste(zOTUbeta$Var1, zOTUbeta$Var2, sep="_")
# Rename columns by replacing "value" with "zOTU"
colnames(zOTUbeta)[3:ncol(zOTUbeta)] <- gsub("value", "zOTU", colnames(zOTUbeta)[3:ncol(zOTUbeta)])
colnames(otu)[3:ncol(otu)] <- gsub("value", "OTU", colnames(otu)[3:ncol(otu)])
zOTUbeta <- merge(zOTUbeta[,3:ncol(zOTUbeta)], otu, by="index")

# Join on size data, because we will only use kipuka > 5000 m^2
zOTUbeta <- merge (zOTUbeta, richness[,c(1,9,10)], by.x="Var1", by.y="ID", all.x=T)
zOTUbeta <- merge (zOTUbeta, richness[,c(1,9,10)], by.x="Var2", by.y="ID", all.x=T)

# We need a category for coloring our box and whisker plots... 
zOTUbeta <- zOTUbeta %>%
  mutate(Site = case_when(
    Site.x == "Lava" & Site.y == "Lava" ~ "Lava",
    Site.x == "Center" & Site.y == "Center" ~ "Center",
    Site.x == "Edge" & Site.y == "Edge" ~ "Edge",
    Site.x == "Stainback" & Site.y == "Stainback" ~ "Stainback",
    Site.x == "Kona" & Site.y == "Kona" ~ "Kona",
    (Site.x == "Center" & Site.y == "Edge") | (Site.x == "Edge" & Site.y == "Center") ~ "Center-Edge",
    (Site.x == "Center" & Site.y == "Stainback") | (Site.x == "Stainback" & Site.y == "Center") ~ "Stainback-Center",
    (Site.x == "Center" & Site.y == "Kona") | (Site.x == "Kona" & Site.y == "Center") ~ "Kona-Center",
    (Site.x == "Center" & Site.y == "Lava") | (Site.x == "Lava" & Site.y == "Center") ~ "Lava-Center",
    (Site.x == "Kona" & Site.y == "Edge") | (Site.x == "Edge" & Site.y == "Kona") ~ "Kona-Edge",
    (Site.x == "Kona" & Site.y == "Stainback") | (Site.x == "Stainback" & Site.y == "Kona") ~ "Kona-Stainback",
    (Site.x == "Kona" & Site.y == "Lava") | (Site.x == "Lava" & Site.y == "Kona") ~ "Kona-Lava",
    (Site.x == "Stainback" & Site.y == "Edge") | (Site.x == "Edge" & Site.y == "Stainback") ~ "Stainback-Edge",
    (Site.x == "Stainback" & Site.y == "Lava") | (Site.x == "Lava" & Site.y == "Stainback") ~ "Stainback-Lava",
    (Site.x == "Lava" & Site.y == "Edge") | (Site.x == "Edge" & Site.y == "Lava") ~ "Lava-Edge",
  ))
zOTUbeta$Site <- factor(zOTUbeta$Site, levels=c("Lava", "Edge", "Center", "Stainback", "Kona", 
                                          "Edge-Center", "Lava-Center", "Stainback-Center", "Kona-Center", "Center-Edge",        
                                          "Lava-Edge", "Stainback-Edge", "Kona-Edge", "Stainback-Lava", "Kona-Lava", "Kona-Stainback"))  
zOTUbeta <- zOTUbeta[!is.na(zOTUbeta$index), ]

zOTUbeta <- zOTUbeta[zOTUbeta$Var2 != zOTUbeta$Var1,]

zOTUbeta <- zOTUbeta %>%
  arrange(Var2, Var1)

#########################################################
# Now, compare distance decay of dissimilarity between continuous forest and kipuka habitat 

# Subset only within-habitat-type comparisons
acari_beta <- zOTUbeta[zOTUbeta$Site.x == zOTUbeta$Site.y,]
# We only want continuous forest and kipuka sites for this analysis... 
acari_beta <- acari_beta[acari_beta$Site.x %in% c("Center", "Edge", "Stainback", "Kona"),]
acari_beta$group[acari_beta$Site.x=="Center" | acari_beta$Site.x=="Edge"] <- "Kipuka"
acari_beta$group[acari_beta$Site.x=="Kona" | acari_beta$Site.x=="Stainback"] <- "Continuous forest"

# Join on the distance data
acari_beta <- merge(acari_beta, geo_dist, by="index")


#################################################################################    
# Now plot it

# First, zOTU :
for (i in (1:6)){
LETTER <- c("A.", "B.", "C.", "D.", "E.", "F.")[i]
letter  <- c("a", "b", "c", "d", "e", "f")[i]
ORDER <- c("Hemiptera_zOTU", "Araneae_zOTU", "Psocoptera_zOTU", "Lepidoptera_zOTU", "Coleoptera_zOTU", "Diptera_zOTU")[i] 
order <- c("Hemiptera", "Araneae", "Psocoptera", "Lepidoptera", "Coleoptera", "Diptera")[i]         
plot <- ggplot(data = acari_beta[acari_beta$group == "Kipuka",], 
                 aes(x = dist, y = !!sym(ORDER), colour = Site.x)) +
  geom_smooth(method = 'lm', formula = y ~ log10(x), se = FALSE, size = 1, alpha = 0.2) +
  geom_point(alpha = 0.70, size = 7, shape = 0) +
  scale_colour_manual(values = SiteColors) +
  labs(title = paste0(LETTER, "        ", order), y = "3% OTU beta diversity") +
  KipukaTheme +
  guides(color = "none", shape = "none", fill = "none", linetype = "none") +
  scale_y_continuous(limits = c(0.3, 1), breaks=seq(0, 1, 0.4)) +
  scale_x_continuous(expand = c(0.1, 0.1), breaks=seq(0, 10000, 2000)) +
  theme(strip.text = element_text(size = 70),
        panel.grid.major = element_line(
          rgb(105, 105, 105, maxColorValue = 255),
          linetype = "dotted",
          size = 1),
        panel.grid.minor = element_line(
          rgb(105, 105, 105, maxColorValue = 255),
          linetype = "dotted",
          size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 50),
        plot.title = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.title = element_text(size = 50),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 15))
assign(letter, plot)        
}
# Save the plot
jpeg("../Figures/Order_zOTU_Jaccard.jpg", width = 2000, height = 2000)
plot_grid(a, b, c, d, e, f, nrow = 3, ncol=2) +  
draw_label("Distance (m)", x = 0.5, y=0, vjust = -0.5, size = 50, fontfamily = "serif") +
draw_label("zOTU beta diversity", x = 0, y = 0.5, vjust = 1, angle = 90, size = 50, fontfamily = "serif") 
dev.off()

# Now, 3% radius OTU :
for (i in (1:6)){
LETTER <- c("A.", "B.", "C.", "D.", "E.", "F.")[i]
letter  <- c("a", "b", "c", "d", "e", "f")[i]
ORDER <- c("Hemiptera_OTU", "Araneae_OTU", "Psocoptera_OTU", "Lepidoptera_OTU", "Coleoptera_OTU", "Diptera_OTU")[i] 
order <- c("Hemiptera", "Araneae", "Psocoptera", "Lepidoptera", "Coleoptera", "Diptera")[i]         
plot <- ggplot(data = acari_beta[acari_beta$group == "Kipuka",], 
                 aes(x = dist, y = !!sym(ORDER), colour = Site.x)) +
  geom_smooth(method = 'lm', formula = y ~ log10(x), se = FALSE, size = 1, alpha = 0.2) +
  geom_point(alpha = 0.70, size = 7, shape = 15) +
  scale_colour_manual(values = SiteColors) +
  labs(title = paste0(LETTER, "        ", order), y = "3% OTU beta diversity") +
  KipukaTheme +
  guides(color = "none", shape = "none", fill = "none", linetype = "none") +
  scale_y_continuous(limits = c(0.1, 1), breaks=seq(0, 1, 0.4)) +
  scale_x_continuous(expand = c(0.1, 0.1), breaks=seq(0, 10000, 2000)) +
  theme(strip.text = element_text(size = 70),
        panel.grid.major = element_line(
          rgb(105, 105, 105, maxColorValue = 255),
          linetype = "dotted",
          size = 1),
        panel.grid.minor = element_line(
          rgb(105, 105, 105, maxColorValue = 255),
          linetype = "dotted",
          size = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 50),
        plot.title = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.title = element_text(size = 50),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 12))
assign(letter, plot)        
}
# Save the plot
jpeg("../Figures/Order_3OTU_BrayCurtis.jpg", width = 2000, height = 2000)
plot_grid(a, b, c, d, e, f, nrow = 3, ncol=2) +  
draw_label("Distance (m)", x = 0.5, y=0, vjust = -0.5, size = 50, fontfamily = "serif") +
draw_label("3% OTU beta diversity", x = 0, y = 0.5, vjust = 1, angle = 90, size = 50, fontfamily = "serif") 
dev.off()

