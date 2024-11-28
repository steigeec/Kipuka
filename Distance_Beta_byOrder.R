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







































richness <- read.csv("merged_by_site_2.csv")
dist_beta <- read.csv("Distance_v_beta.csv")
dist_diff <- read.csv("Distance_v_differentiation.csv")
geo_dist<-read.csv("geo_dist.csv")
otu <- read.csv("OTUs.csv")

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

#Create a function for plotting R2
eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                list(a = format(coef(m)[1], digits = 4),
                b = format(coef(m)[2], digits = 4),
                r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}
                   
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
#geo_dist$log_dist<-log(geo_dist$dist+0.00001)
geo_dist<-geo_dist[,c(3, 4)]
names(geo_dist)<-c("geo_dist", "index")
#xy <- t(combn(colnames(geo_dist), 2))
#geo_dist <- data.frame(xy, dist=geo_dist[xy])

OTU <- otu[8:67, 29:3051] 
#Make first column the row names
OTU["13",1]<-"ZotuID"
OTU["14",1]<-"threepOTU"
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)

#What are the orders for which we need to tailor?
unique(OTU$Order)
#Classes: Acari
#Order: (1) "Coleoptera"       (2) "Diptera"          (3) "Hemiptera"        (4) "Lepidoptera"        (5) "Psocoptera"      (6) Aranea 

orders<-c("Diptera", "Hemiptera", "Lepidoptera", "Psocoptera", "Araneae", "Coleoptera")
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
        acari_beta <- vegdist(acari, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
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
order_all<-rbind(Araneae_beta, Diptera_beta, Hemiptera_beta, Lepidoptera_beta, Psocoptera_beta, Coleoptera_beta) 

#########################################################################################
# PLOT BETA DIVERSITY BY DISTANCE FOR EACH ORDER
# THESE PLOTS ARE NOT IN THE PAPER

jpeg("../Figures/Order_beta_diversity.jpg", width=1500, height=2000)
ggplot(data=order_all) + 
  geom_smooth(method='lm', aes(x=geo_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20, shape=15)+
  geom_point(aes(x=geo_dist, y=dist, colour=Site.x), alpha=0.70, size=4, shape=18) + 
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual("Site type", values=SiteColors) +
  labs(title="", x="Distance (m)", y="zOTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=3)+
  KipukaTheme +
  coord_cartesian(ylim=c(0.4, 1))+
  guides(colour="none")+
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +                   
  #geom_text(x = min(order_all$geo_dist), y = max(order_all$dist), label = eq(order_all$geo_dist,order_all$dist), parse = TRUE)+                                      
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title.x=element_text(size=45, margin=margin(-15,0,0,0)), 
       axis.title.y=element_text(size=45, margin=margin(0,-15,0,0)), 
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()

#Save another version of this jpeg, now only highlighting edge and center turnover
jpeg("../Figures/Order_beta_diversity_kipukas_zOTU.jpg", width=1500, height=2000)
ggplot(data=order_all[order_all$Site.x==c("Center", "Edge"),]) + 
  geom_smooth(method='lm', aes(x=geo_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20)+
  geom_point(aes(x=geo_dist, y=dist, colour=Site.x), alpha=0.70, size=4, shape=0, stroke=3.5) + 
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual("Position in kipuka", values=SiteColors, limits = c("Center", "Edge")) +
  labs(title="", x="Log distance (m)", y="zOTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=3)+
  KipukaTheme +
  coord_cartesian(ylim=c(0.4, 1))+
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
       axis.title.x=element_text(size=45, margin=margin(-15,0,0,0)), 
       axis.title.y=element_text(size=45, margin=margin(0,-15,0,0)), 
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()

#######################################################################################
# NOW, PLOT AGAIN-- THIS TIME REMOVING THE 5 SMALLEST KIPUKA THAT BEHAVE JUST LIKE EDGE
                     
to_rem <- c("1K01C", "1K01E", "1K06C", "1K06E", "1K07C", "1K07E", "1K12C", "1K12E", "1K13C", "1K13E")
sub_all <- order_all[!(order_all$row %in% to_rem) & !(order_all$col %in% to_rem),]                     

for (i in 1:length(unique(sub_all$order))){  
    for (j in 1:length(c("Edge", "Center"))){
        taxon<-unique(sub_all$order)[i]
        area<-c("Edge", "Center")[j]
        print(paste0("Results for ",taxon," and ",area))                   
        gam_model <- gam(dist ~ s(log10(geo_dist)), data = sub_all[sub_all$order==taxon,])
        # Check assumptions
        # 1. Residuals vs Fitted Values Plot
        par(mar = c(1, 1, 1, 1))                         
        plot(residuals(gam_model) ~ fitted(gam_model), main = "Residuals vs Fitted", xlab = "Fitted Values",ylab = "Residuals")
        abline(h = 0, col = "red", lty = 2)
        # 2. Normal Q-Q Plot
        qqnorm(residuals(gam_model))
        qqline(residuals(gam_model), col = "red")
        # 3. Scale-Location (Spread-Location) Plot
        plot(sqrt(abs(residuals(gam_model))) ~ fitted(gam_model), main = "Scale-Location Plot", xlab = "Fitted Values", ylab = "sqrt(|Residuals|)")
        abline(h = 0, col = "red", lty = 2)
        # Print summary of the linear model
        print(summary(gam_model))    
}
}
                     
jpeg("../Figures/Order_beta_diversity_kipukas_sub_zOTU.jpg", width=1500, height=2000)
ggplot(data=sub_all[sub_all$Site.x==c("Center", "Edge"),], aes(x=geo_dist, y=dist, colour=Site.x)) + 
  geom_smooth(method = "lm", formula = y ~ log10(x), se = F, size=1, alpha=0.20, aes(fill=Site.x)) +                   
  geom_point(alpha=0.70, size=4, shape=0, stroke=2) + 
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual("Position in kipuka", values=SiteColors, limits = c("Center", "Edge")) +
  labs(title="", x="Distance (m)", y="zOTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=3)+
  KipukaTheme +
  #coord_cartesian(ylim=c(0.4, 1))+
  guides(colour="none")+
  #scale_x_continuous(trans='log10',
  #                   breaks=trans_breaks('log10', function(x) 10^x),
  #                   labels=trans_format('log10', math_format(10^.x)))  +
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
       legend.position = "top",
       plot.margin = margin(1, 35, 1, 1))  
dev.off()                     
                     
#################################################################################
#Repeat the same as above, now for method jaccard

orders<-c("Diptera", "Hemiptera", "Lepidoptera", "Psocoptera", "Araneae", "Coleoptera")
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
        acari_beta <- vegdist(acari, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
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
order_all<-rbind(Araneae_beta, Coleoptera_beta, Diptera_beta, Hemiptera_beta, Lepidoptera_beta, Psocoptera_beta, Coleoptera_beta) 
to_rem <- c("1K01C", "1K01E", "1K06C", "1K06E", "1K07C", "1K07E", "1K12C", "1K12E", "1K13C", "1K13E")
sub_all <- order_all[!(order_all$row %in% to_rem) & !(order_all$col %in% to_rem),]                     

##############
# First, let's get the stats for these relationships:
# First, test assumptions:  
for (i in 1:length(unique(sub_all$order))){  
    for (j in 1:length(c("Edge", "Center"))){
        taxon<-unique(sub_all$order)[i]
        area<-c("Edge", "Center")[j]
        print(paste0("Results for ",taxon," and ",area))                   
        gam_model <- gam(dist ~ s(log10(geo_dist)), data = sub_all[sub_all$order==taxon,])
        # Check assumptions
        # 1. Residuals vs Fitted Values Plot
        par(mar = c(1, 1, 1, 1))                         
        plot(residuals(gam_model) ~ fitted(gam_model), main = "Residuals vs Fitted", xlab = "Fitted Values",ylab = "Residuals")
        abline(h = 0, col = "red", lty = 2)
        # 2. Normal Q-Q Plot
        qqnorm(residuals(gam_model))
        qqline(residuals(gam_model), col = "red")
        # 3. Scale-Location (Spread-Location) Plot
        plot(sqrt(abs(residuals(gam_model))) ~ fitted(gam_model), main = "Scale-Location Plot", xlab = "Fitted Values", ylab = "sqrt(|Residuals|)")
        abline(h = 0, col = "red", lty = 2)
        # Print summary of the linear model
        print(summary(gam_model))    
}
}
        
jpeg("../Figures/Order_beta_diversity_JACCARD.jpg", width=1500, height=2000)
ggplot(data=sub_all[sub_all$Site.x==c("Center", "Edge"),], aes(x=geo_dist, y=dist, colour=Site.x)) + 
  geom_smooth(method = "lm", formula = y ~ log10(x), se = F, size=1, alpha=0.20, aes(fill=Site.x)) +                   
  geom_point(alpha=0.70, size=4, shape=0, stroke=2) + 
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual("Site type", values=SiteColors, limits = c("Center", "Edge")) +
  labs(title="", x="Distance (m)", y="3 % radius OTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=4)+
  KipukaTheme +
  #coord_cartesian(ylim=c(0.6, 1))+
  guides(colour="none")+
  #scale_x_continuous(trans='log10',
  #                   breaks=trans_breaks('log10', function(x) 10^x),
  #                   labels=trans_format('log10', math_format(10^.x)))  +
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
       legend.position = "top", 
        plot.margin = margin(1, 35, 1, 1))  
dev.off()
                  

to_rem <- c("1K01C", "1K01E", "1K06C", "1K06E", "1K07C", "1K07E", "1K12C", "1K12E", "1K13C", "1K13E")
sub_all <- order_all[!(order_all$row %in% to_rem) & !(order_all$col %in% to_rem),] 
                     
#Save an alternate version of this jpeg, highlighting only kipuka data
jpeg("../Figures/Order_beta_diversity_JACCARD_kipukas_zOTU.jpg", width=1500, height=2000)
ggplot(data=sub_all[sub_all$Site.x==c("Center", "Edge"),]) + 
  geom_smooth(method='lm', aes(x=geo_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20)+
  geom_point(aes(x=geo_dist, y=dist, colour=Site.x), alpha=0.70, size=4, stroke=2.5, shape=0) + 
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual("Position in kipuka", values=SiteColors, limits = c("Center", "Edge")) +
  labs(title="", x="Log distance (m)", y="zOTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=4)+
  KipukaTheme +
  coord_cartesian(ylim=c(0.6, 1))+
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
       axis.title.x=element_text(size=45, margin=margin(-15,0,0,0)), 
       axis.title.y=element_text(size=45, margin=margin(0,-15,0,0)), 
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()
##################################################################################################################################
#NOW, FOR 3% OTU BETA DIVERSITY

##HERE IS WHERE I'm at... NEED TO POOL FOR 3%OTU
OTU[] <- lapply(OTU, function(x) if(is.character(x)){
              factor(trimws(x))
              } else x
        )  
OTU$threepOTU<-as.factor(OTU$threepOTU)
#Can't be na in 3%otu col
OTU<-OTU[!(is.na(OTU$threepOTU) | OTU$threepOTU==""), ]
OTU[9:60]<-sapply(OTU[9:60],as.numeric)     
                   
otu3<-OTU[,c(7,9:60)]%>%
        group_by(threepOTU)%>%
        mutate(sum = rowSums(across(where(is.numeric)), na.rm=TRUE))
        #summarize(Class="Class", Order="Order", across(9:60, sum))

for (ORDER in 1:length(orders)){
        O<-orders[ORDER]
        acari <- OTU[OTU$Order==O,] 
        #get rid of empty rows
        acari <- acari[rowSums(is.na(acari)) != ncol(acari), ]  
        acari<-acari[,9:60]
        #Sites must be rows, and species are columns
        acari<-as.data.frame(t(as.matrix(acari)))
        acari[] <- lapply(acari, as.numeric)
        acari_beta <- vegdist(acari, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
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
order_all<-rbind(Araneae_beta, Diptera_beta, Hemiptera_beta, Lepidoptera_beta, Psocoptera_beta, Coleoptera_beta) 
to_rem <- c("1K01C", "1K01E", "1K06C", "1K06E", "1K07C", "1K07E", "1K12C", "1K12E", "1K13C", "1K13E")
sub_all <- order_all[!(order_all$row %in% to_rem) & !(order_all$col %in% to_rem),] 

#Save another version of this jpeg, now only highlighting edge and center turnover
jpeg("../Figures/Order_beta_diversity_kipukas_3perc.jpg", width=1500, height=2000)
ggplot(data=sub_all[sub_all$Site.x==c("Center", "Edge"),]) + 
  geom_smooth(method='lm', aes(x=geo_dist, y=dist, colour=Site.x, fill=Site.x), size=1, alpha=0.20)+
  geom_point(aes(x=geo_dist, y=dist, colour=Site.x), alpha=0.70, size=4, shape=15) + 
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge")) +
  scale_fill_manual("Position in kipuka", values=SiteColors, limits = c("Center", "Edge")) +
  labs(title="", x="Distance (m)", y="3% OTU beta diversity") +
  facet_wrap(~order, ncol=2, nrow=4)+
  KipukaTheme +
  coord_cartesian(ylim=c(0.4, 1))+
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
       axis.title.x=element_text(size=45, margin=margin(-15,0,0,0)), 
       axis.title.y=element_text(size=45, margin=margin(0,-15,0,0)), 
        axis.text = element_text(size=40), 
        plot.title=element_text(size=45), 
        legend.text=element_text(size=40), 
        legend.title = element_text(size=40),
       legend.position = "top")
dev.off()
