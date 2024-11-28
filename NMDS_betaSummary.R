
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
ExtendedSiteColors <- c("Edge-Center"="grey24", "Lava-Center"="grey24", "Stainback-Center"="grey24", "Kona-Center"="grey24", "Center-Edge"="grey24", 
                        "Lava-Edge"="grey24", "Stainback-Edge"="grey24", "Kona-Edge"="grey24", "Center"="#332288", "Edge"="#6699CC", 
                        "Lava"="#888888", "Stainback-Lava"="grey24", "Kona-Lava"="grey24", "Stainback"="#999933", "Kona-Stainback"="grey24", "Kona"="#117733")
  
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
#Version with all the data

richness <- read.csv("merged_by_site_2.csv")
names(richness) <- richness[18,]
richness <- richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))
richness$Site<-gsub("Stainbeck","Stainback",as.character(richness$Site))
# beta diversity for 3 % OTU for within-area
otu <- read.csv("OTU3_Bray.csv")
# beta diversity for zOTU for within-area
zOTUbeta <- read.csv("zOTU_jaccard.csv")
# geographic distances
geo_dist<-read.csv("geo_dist.csv")
nmds<-read.csv("nmds3otu.csv")
nmds$Site<-gsub("Stainbeck","Stainback",as.character(nmds$Site))

######################################################33

#Wrange geo distances
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

#######################

# merge zOTU and 3OTU data
otu$index <- paste(otu$Var1, otu$Var2, sep="_")
zOTUbeta$index <- paste(zOTUbeta$Var1, zOTUbeta$Var2, sep="_")
colnames(zOTUbeta)[3] <- "zOTU"
colnames(otu)[3] <- "OTU"
zOTUbeta <- merge(zOTUbeta[,3:4], otu, by="index")

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

##########################
# With the table of within-area-type dissimilarity measures,  compare the dissimilarity within small edges to that within big edges,
# and the dissimilarity of small cores to that within big cores...

# Remove anything that has something other than Center or Edge in Site.x or Site.y
size <- zOTUbeta[zOTUbeta$Site.x %in% c("Edge") & zOTUbeta$Site.y %in% c("Edge") | zOTUbeta$Site.x %in% c("Center") & zOTUbeta$Site.y %in% c("Center"), ]
size$Area.x<-as.numeric(size$Area.x)
size$Area.y<-as.numeric(size$Area.y)

length(unique(size$Area.x[size$Area.x > 5000]))   # 8 large
length(unique(size$Area.x[size$Area.x < 5000]))   # 5 small

# subset 8 largest and the 4 smallest kipuka
size$group[size$Area.x > 5000 & size$Area.y > 5000]<-"large"
size$group[size$Area.x < 5000 & size$Area.y < 5000]<-"small"

# Remove flipped pairs... 
size <- size %>%
  distinct(OTU, .keep_all = TRUE)

# Conduct a formal test.   Do "large-large" have greater beta than "small-small" ?

###### 3 % radius OTU
# Step 1: Test for Assumptions
# Shapiro-Wilk test for normality
shapiro.test(size$OTU[size$group == "large" & size$Site.x=="Center"])   # p-value = 0.3497
shapiro.test(size$OTU[size$group == "small" & size$Site.x=="Center"])   # p-value = 0.3401
# Levene's test for homogeneity of variances   
leveneTest(size$OTU[size$group %in% c("large", "small")& size$Site.x=="Center"] ~ size$group[size$group %in% c("large", "small")& size$Site.x=="Center"])     # 0.7383
# Step 2: Conduct the Parametric Test (t-test)
# Assuming assumptions are met
# Perform independent t-test
wilcox.test(size$OTU[size$group %in% c("large", "small") & size$Site.x == "Center"] ~ 
            size$group[size$group %in% c("large", "small") & size$Site.x == "Center"])
shapiro.test(size$OTU[size$group == "large" & size$Site.x=="Edge"])   # p-value = 0.3663
shapiro.test(size$OTU[size$group == "small" & size$Site.x=="Edge"])   # p-value = 0.6378
# Levene's test for homogeneity of variances   
leveneTest(size$OTU[size$group %in% c("large", "small")& size$Site.x=="Edge"] ~ size$group[size$group %in% c("large", "small")& size$Site.x=="Edge"])     # 0.1295
# Step 2: Conduct the Parametric Test (t-test)
# Assuming assumptions are met
# Perform independent t-test
wilcox.test(size$OTU[size$group %in% c("large", "small") & size$Site.x == "Edge"] ~ 
            size$group[size$group %in% c("large", "small") & size$Site.x == "Edge"])

###### zOTU
# Step 1: Test for Assumptions
# Shapiro-Wilk test for normality
shapiro.test(size$zOTU[size$group == "large" & size$Site.x=="Center"])   # p-value = 0.3497
shapiro.test(size$zOTU[size$group == "small" & size$Site.x=="Center"])   # p-value = 0.3401
# Levene's test for homogeneity of variances   
leveneTest(size$zOTU[size$group %in% c("large", "small")& size$Site.x=="Center"] ~ size$group[size$group %in% c("large", "small")& size$Site.x=="Center"])     # 0.7383
# Step 2: Conduct the Parametric Test (t-test)
# Assuming assumptions are met
# Perform independent t-test
wilcox.test(size$zOTU[size$group %in% c("large", "small") & size$Site.x == "Center"] ~ 
            size$group[size$group %in% c("large", "small") & size$Site.x == "Center"])
shapiro.test(size$zOTU[size$group == "large" & size$Site.x=="Edge"])   # p-value = 0.3663
shapiro.test(size$zOTU[size$group == "small" & size$Site.x=="Edge"])   # p-value = 0.6378
# Levene's test for homogeneity of variances   
leveneTest(size$zOTU[size$group %in% c("large", "small")& size$Site.x=="Edge"] ~ size$group[size$group %in% c("large", "small")& size$Site.x=="Edge"])     # 0.1295
# Step 2: Conduct the Parametric Test (t-test)
# Assuming assumptions are met
# Perform independent t-test
wilcox.test(size$zOTU[size$group %in% c("large", "small") & size$Site.x == "Edge"] ~ 
            size$group[size$group %in% c("large", "small") & size$Site.x == "Edge"])

#######################################################
# Test difference between area size and zOTU and 3% radius OTU turnover

# from zOTUbeta, keep only the comparisons within a single kipuka...
size <- zOTUbeta[zOTUbeta$Site.x %in% c("Edge") & zOTUbeta$Site.y %in% c("Center") | zOTUbeta$Site.x %in% c("Center") & zOTUbeta$Site.y %in% c("Edge"), ]
size <- size[size$Area.x == size$Area.y,]

# Remove flipped pairs... 
size <- size %>%
  distinct(OTU, .keep_all = TRUE)

# First, test assumptions:  
# Fit linear regression model
gam_model <- gam(OTU ~ s(log10(Area.x)), data = size)
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
# 4. Residuals vs Leverage Plot (Cook's distance)
plot(hatvalues(gam_model), cooks.distance(gam_model), main = "Residuals vs Leverage", xlab = "Leverage", ylab = "Cook's distance")
abline(h = 4/length(CE[CE$metric=="3% OTU",]$value), col = "red", lty = 2)
# Print summary of the linear model
summary(gam_model)                                    
                         
######################################################                         
# Shown as a log10 function
jpeg("../Figures/turnover-edge-center_v3.jpg", width = 1000, height = 1000)  # Increase the width                       
ggplot(data = size, aes(x = Area.x, y = OTU)) + 
  geom_smooth(method = "lm", formula = y ~ poly(log10(x), 2), se = FALSE, size=1, alpha=0.20, col="black") +
  geom_point(colour = "#6699CC", fill = "#332288", alpha = 0.70, size = 8, shape = 21, stroke = 7) + 
  labs(x = "Kipuka area (" ~ m^2 ~ ")", y = "3% OTU beta diversity") +
  KipukaTheme +
  # coord_cartesian(ylim = c(0.4, 1)) +
  guides(colour = "none") +
  scale_x_continuous(trans = 'identity', expand = c(0.05, 0.02))  +  # Adjust the expand parameter
  scale_y_continuous(expand = c(0.005, 0.05))  +  # Adjust the expand parameter
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
          rgb(105, 105, 105, maxColorValue = 255),
          linetype = "dotted", 
          size = 1),   
        panel.grid.minor = element_line(
          rgb(105, 105, 105, maxColorValue = 255),
          linetype = "dotted", 
          size = 0.5), 
        axis.title.y = element_text(size = 45), 
        axis.title.x = element_text(size = 45, margin = margin(-25, 0, 0, 0)), 
        axis.text = element_text(size = 40), 
        plot.title = element_text(size = 45), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(size = 40),
        legend.position = "top",
        plot.margin = margin(1, 15, 1, 1))  
dev.off()


#################################################################################

# PERMANOVA 
# to look for overall community comp differences between sites

# given the differences we have demonstrated for small (<5000m2) and large kipuka, 
# at this point remove the smallest kipuka
zOTUbeta <- zOTUbeta[
  (is.na(zOTUbeta$Area.x) & is.na(zOTUbeta$Area.y)) | 
  (zOTUbeta$Area.x > 5000 & zOTUbeta$Area.y > 5000) | 
    (zOTUbeta$Area.x > 5000 & is.na(zOTUbeta$Area.y)) | 
    (is.na(zOTUbeta$Area.x) & zOTUbeta$Area.y > 5000), 
]

IDS <- unique(zOTUbeta$Var1)     
Site <- as.factor(zOTUbeta$Site.x[zOTUbeta$Var1 %in% IDS & zOTUbeta$Var2=="1K02C"])

# FIRST, 3% RADIUS OTU... ###        
# Transform zOTUbeta back into distance matrices... 
distance_matrix_df <- zOTUbeta[c(1,2,5)] %>%
  pivot_wider(names_from = Var2, values_from = OTU) 

# Set rownames from Var1, and then convert to a matrix
distance_matrix_df <- distance_matrix_df %>%
  column_to_rownames(var = "Var1") %>%
  as.matrix()

adonis2(distance_matrix_df ~ Site, permutations = 999)              
# Now, the pairwise test:                 
pairwise.adonis(distance_matrix_df, Site, p.adjust.m="fdr")
                         
# NEXT, ZOTU TESTS... #############################################
# Transform zOTUbeta back into distance matrices... 
distance_matrix_df <- zOTUbeta[c(1,2,4)] %>%
  pivot_wider(names_from = Var2, values_from = zOTU) 

# Set rownames from Var1, and then convert to a matrix
distance_matrix_df <- distance_matrix_df %>%
  column_to_rownames(var = "Var1") %>%
  as.matrix()

adonis2(distance_matrix_df ~ Site, permutations = 999)                        
# Now the pairwise test                       
pairwise.adonis(distance_matrix_df, Site, p.adjust.m="fdr")


#######################################################
# Plot these

zOTUbeta <- zOTUbeta[zOTUbeta$Var1 != zOTUbeta$Var2, ]

#NMDS plot for 3%OTU
#Create an NMDS plot with columns MDS1 and MDS2
nmds$Area<-as.numeric(gsub(",","",as.character(nmds$Area)))
nmds$pointsize<-round(sqrt(as.numeric(nmds$Area))/10,0)
nmds$pointsize[nmds$Site=="Lava" & is.na(nmds$pointsize)] <- 4
nmds$pointsize[nmds$Site=="Kona" & is.na(nmds$pointsize)] <- 4
nmds$pointsize[nmds$Site=="Stainback" & is.na(nmds$pointsize)] <- 4
nmds <- nmds[nmds$Site != "1K08E" & nmds$Site != "1K08C",]                          
                         
a <- ggplot() + 
  geom_point(data=nmds[nmds$Site!=c("C-E", "C-F"),],aes(x=MDS1OTU,y=MDS2OTU,colour=Site, size=pointsize, shape=Site), alpha=0.70, stroke=3) + 
  # Confidence ellipses
  stat_ellipse(
    data = nmds[nmds$Site != c("C-E", "C-F"), ], aes(x = MDS1OTU, y = MDS2OTU, colour = Site), level = 0.95,  # 95% confidence interval
    linetype = "dashed", size=1.5) +
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge", "Lava", "Kona", "Stainback")) +
  scale_shape_manual("Site", values=c("Center" = 16, "Edge" = 16, "Lava"=3, "Kona"=2, "Stainback"=2)) +
  scale_size_continuous("Kipuka area ("~m^2~")", range=c(2,32), breaks=seq(2,32,5), labels=round((10*seq(2,32,5))^2,100)) +
  labs(title="A.", x="NMDS1", y="NMDS2") +
  #coord_equal() +
  scale_x_continuous(breaks=seq(-2,1.5,0.5), limits=c(-1.75, 1.3)) +
  scale_y_continuous(limits=c(-.9, 1)) +
KipukaTheme +
  theme(axis.title.x = element_text(size=80), 
        axis.text.x = element_text(angle=45, size=80, hjust=1, vjust=1),
        legend.position = "none", 
        panel.grid.major = element_line(rgb(105, 105, 105, maxColorValue = 255), linetype = "dotted", size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5),
    plot.title = element_text(size = 70, hjust=-.1), # Offset title to the left
        #plot.margin = margin(0, 0, 1.5, 0, "cm")
       ) # Increased bottom margin

                         
b<- ggplot() + 
  geom_boxplot(data=zOTUbeta,aes(x=Site, y=OTU, fill=Site), color="black", size=1, outlier.size=2.5)+
  #facet_wrap(~metric, scales="free") +
  scale_fill_manual(values=ExtendedSiteColors) +
  labs(title="B.", y="3% OTU beta diversity", x="") +
  KipukaTheme +
  guides(fill="none")+                       
  theme(axis.title.x=element_blank(), 
        strip.text = element_text(size = 80), 
        axis.text.x = element_text(angle=45, size=80, hjust=1, vjust=1),
        axis.text.y = element_text(angle=45, size=80, margin=margin(0,-10,0,0)), 
        axis.title.y = element_text(size = 80, margin=margin(0,-25,0,0)),
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
    plot.title = element_text(size = 80, hjust=-.1), # Offset title to the left
        legend.text=element_text(size=60), 
        legend.title = element_blank(),
       legend.position = "none", 
    #plot.margin = margin(0.2, 1, 1.5, .8, "cm")
       ) # Increased bottom margin
                         
reserved <- ggplot() + 
  geom_point(data=nmds[nmds$Site!=c("C-E", "C-F"),],aes(x=MDS1zOTU,y=MDS2zOTU,colour=Site, size=pointsize, shape=Site), alpha=0.70, stroke=3) + 
  # Confidence ellipses
  stat_ellipse(
    data = nmds[nmds$Site != c("C-E", "C-F"), ], aes(x = MDS1zOTU, y = MDS2zOTU, colour = Site),
    level = 0.95,  # 95% confidence interval
    linetype = "dashed", size=1.5) +
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge", "Lava", "Kona", "Stainback")) +
  scale_shape_manual("Site", values=c("Center" = 16, "Edge" = 16, "Lava"=3, "Kona"=2, "Stainback"=2)) +
  scale_size_continuous("Kipuka area ("~m^2~")", range=c(2,32), breaks=seq(2,32,5), labels=round((10*seq(2,32,5))^2,100)) +
  labs(title="C.", x="NMDS1", y="NMDS2") +
  #coord_equal() +
  scale_x_continuous(breaks=seq(-2,1.5,0.5), limits=c(-1.5, 1.5)) +
  scale_y_continuous(limits=c(-1, 1)) +
guides(
    fill = guide_legend(override.aes = list(shape = c(21, 22, 23, 24, 25), size = 4)),
    shape = "none"
  ) + 
KipukaTheme +
  theme(       legend.position = "none", 
        legend.title = element_text(size=50),
        legend.text=element_text(size=45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5),
        plot.title=element_text(size=80, hjust=-.1))                     
                     
d<- ggplot() + 
  geom_boxplot(data=zOTUbeta,aes(x=Site, y=zOTU, fill=Site), color="black", size=1, outlier.size=2.5)+
  scale_fill_manual(values=ExtendedSiteColors) +
  labs(title="C.", y="zOTU beta diversity", x="") +
  KipukaTheme +
  guides(fill="none")+                       
  theme(strip.text = element_text(size = 80), 
        axis.text.x = element_text(angle=45, size=80, hjust=1, vjust=1), 
        axis.text.y = element_text(angle=45, size=80, margin=margin(0,-10,0,0)), 
        axis.title.y = element_text(size = 80, margin=margin(0,-25,0,0)),
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=80), 
        plot.title=element_text(size=80, hjust=-.1), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top", 
        #plot.margin = margin(0.2,1,0,.8, "cm")
       )  
                         
jpeg("../Figures/NMDS-turnovers.jpg", width=5000, height=1500) 
plot_grid(a,b,d, nrow=1, ncol=3, rel_widths=c(1.75, 2, 2))                         
dev.off()                                 

#######################################################
# Test for differences between
for (i in 1:2) {
    level <- c("OTU", "zOTU")[i]  # Get the column name as a string
    print(level)
    
    # Dynamically subset the data for the current level
    zOTUbeta$level_value <- zOTUbeta[[level]]  # Create a temporary column
    
    # First, test assumptions
    residuals <- lm(level_value ~ Site, data = zOTUbeta)$residuals
    
    # Plot Normal Q-Q Plot
    par(mar = c(1, 1, 1, 1))
    qqPlot(residuals, main = paste("Normal Q-Q Plot of Residuals for", level))
    
    # Check homogeneity of variances
    p1 <- leveneTest(level_value ~ Site, data = zOTUbeta)[1, 3]
    print(paste0("p-value for Levene's is ", p1))
    
    # Shapiro-Wilk test for normality (optional)
    p2 <- shapiro.test(residuals)$p.value
    print(paste0("p-value for Shapiro-Wilk's is ", p2))
    
    if (p1 < 0.05 || p2 < 0.05) {   
        # If assumptions are violated, conduct Kruskal-Wallis test
        kruskal_result <- kruskal.test(level_value ~ Site, data = zOTUbeta)
        cat("Kruskal-Wallis test results for", level, "\n")
        print(kruskal_result)
        
        # Pairwise Wilcoxon test
        pairwise_result <- pairwise.wilcox.test(zOTUbeta$level_value, zOTUbeta$Site, p.adj = "bonferroni")
        cat("Pairwise Wilcoxon test results for", level, "\n")
        print(pairwise_result)
    } else {   
        # Conduct ANOVA
        anova_result <- aov(level_value ~ Site, data = zOTUbeta)
        cat("ANOVA test results for", level, "\n")
        print(summary(anova_result))
        
        # Post hoc Tukey's HSD test
        tukey_result <- TukeyHSD(anova_result)
        cat("Tukey's HSD test results for", level, "\n")
        print(tukey_result)
    }
}

# Clean up the temporary column
zOTUbeta$level_value <- NULL



















##########################
#Join together tables of each part...?
dist_beta <- dist_beta_0[,1:7]
i<-c(3, 4, 5, 6, 7)
dist_beta[ , i] <- apply(dist_beta[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))

Center <- dist_beta[!is.na(dist_beta$Center),]
Center$site <- "Center"
names(Center)[names(Center) == "Center"] <- "beta"
Center <- Center[ , colSums(is.na(Center)) < nrow(Center)]                    
                         
Stainbeck <- dist_beta[!is.na(dist_beta$Stainbeck),]
Stainbeck$site <- "Stainbeck"  
names(Stainbeck)[names(Stainbeck) == "Stainbeck"] <- "beta"
Stainbeck <- Stainbeck[ , colSums(is.na(Stainbeck)) < nrow(Stainbeck)]   
                         
Kona <- dist_beta[!is.na(dist_beta$Kona),]
Kona$site <- "Kona"
names(Kona)[names(Kona) == "Kona"] <- "beta"
Kona <- Kona[ , colSums(is.na(Kona)) < nrow(Kona)]   
                         
Edge <- dist_beta[!is.na(dist_beta$Edge),]
Edge$site <- "Edge"
names(Edge)[names(Edge) == "Edge"] <- "beta"
Edge <- Edge[ , colSums(is.na(Edge)) < nrow(Edge)] 
                         
Lava <- dist_beta[!is.na(dist_beta$Lava),]
Lava$site <- "Lava"
names(Lava)[names(Lava) == "Lava"] <- "beta"
Lava <- Lava[ , colSums(is.na(Lava)) < nrow(Lava)]                          
                         
dist_beta <- rbind(Center, Stainbeck, Kona, Edge, Lava)
names(dist_beta)[names(dist_beta) == "ï..dist"] <- "dist"                       

#Now summarize beta diversity between site types
OTU3beta<-OTU3beta[,c(7, 4, 5)]     # logdist, beta, Site.x
OTU3beta$metric <- "3% OTU"
dist_beta<-dist_beta[,c(2, 3, 4)]       # logdist, beta, site
dist_beta$metric <- "zOTU"
names(dist_beta)<-c("log_dist", "beta", "Site.x", "metric")

beta <- rbind(dist_beta, OTU3beta) 
#Reorder facets
beta$metric <- factor(beta$metric, levels = rev(c("zOTU", "3% OTU")))   
beta <- beta[beta$Site != "1K08E" & beta$Site != "1K08C",]
beta <- beta[beta$Site.x != "C-F" & beta$Site.x != "C-E",]

#Fix spelling error on sheet before proceeding
beta$Site.x<-gsub("Stainbeck","Stainback",as.character(beta$Site.x))
                         
beta$Site.x<-as.factor(beta$Site.x)
beta$Site.x <- factor(beta$Site.x, levels=c("Lava", "Edge", "Center", "Stainback", "Kona"))  

#What is the average zOTU for each site type? 3%OTU?
for (LEVEL in 1:length(levels(beta$Site.x))){
        type<-levels(beta$Site.x)[LEVEL]
        average_zotu<-mean(beta$beta[beta$Site.x==type & beta$metric=="zOTU"])
        average_3otu<-mean(beta$beta[beta$Site.x==type & beta$metric=="3% OTU"])
        print(paste0("The zOTU beta of ",type," is ",average_zotu))
        print(paste0("The 3% OTU beta of ",type," is ",average_3otu))
        }              

#Bray curtis center-edge pairs only
CE<-merge(CE, richness, by.x="log_dist", by.y="ï..ID")

#Fix spelling error on sheet before proceeding
CE$Site<-gsub("Stainbeck","Stainback",as.character(CE$Site))                           
                     
