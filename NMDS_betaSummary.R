#Project: Kipukas               #
#Script: NMDS.R                 #
#Author: Emma Steigerwald       #
#Date:17 Feb 2022               #
#################################

#Set up my working environment
setwd("G:/My Drive/Kipuka/Data")
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

nmds$pointsize<-round(sqrt(as.numeric(nmds$Area))/10,0)
nmds$pointsize[nmds$Site=="Lava" & is.na(nmds$pointsize)] <- 2
nmds$pointsize[nmds$Site=="Kona" & is.na(nmds$pointsize)] <- 2
nmds$pointsize[nmds$Site=="Stainback" & is.na(nmds$pointsize)] <- 2

nmds <- nmds[nmds$Site != "1K08E" & nmds$Site != "1K08C",] 

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
                out_row<-out_row+1}}

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
                out_row<-out_row+1} }

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

#######################
# Produce a table of within-area-type dissimilarity measures:

#Grab 3% OTU beta diversity metrics
acari_beta<-as.matrix(otu)
#Remove all the XXXX from colnames
colnames(acari_beta)<- gsub('[X]', '', colnames(acari_beta))
rownames(acari_beta)<-acari_beta[,1]
acari_beta<-acari_beta[,-1]
acari_beta<-data.frame(col=colnames(acari_beta)[col(acari_beta)], row=rownames(acari_beta)[row(acari_beta)], beta=c(acari_beta))
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
acari_beta <- acari_beta%>% distinct(beta, .keep_all= TRUE)
acari_beta$beta<-as.numeric(acari_beta$beta)

##########################
# With the table of within-area-type dissimilarity measures,  compare the dissimilarity within small edges to that within big edges,
# and the dissimilarity of small cores to that within big cores...

# Join the size data onto acari_beta
size<-merge(acari_beta, richness[c(1,10)], by.x="row", by.y="ï..ID")
size<-merge(size, richness[c(1,10)], by.x="col", by.y="ï..ID")

# Remove anything that has something other than Center or Edge in Site.x or Site.y
size <- size[size$Site.x %in% c("Center", "Edge") & size$Site.y %in% c("Center", "Edge"), ]
size$Area.x<-as.numeric(size$Area.x)
size$Area.y<-as.numeric(size$Area.y)

length(unique(size$Area.x[size$Area.x > 5000]))   # 8 large
length(unique(size$Area.x[size$Area.x < 5000]))   # 4 small

# subset 8 largest and the 4 smallest kipuka
size$group[size$Area.x > 5000 & size$Area.y > 5000]<-"large"
size$group[size$Area.x < 5000 & size$Area.y < 5000]<-"small"

# Conduct a formal test.   Do "large-large" have greater beta than "small-small" ?
# Step 1: Test for Assumptions
# Shapiro-Wilk test for normality
shapiro.test(size$beta[size$group == "large" & size$Site.x=="Center"])   # p-value = 0.3497
shapiro.test(size$beta[size$group == "small" & size$Site.x=="Center"])   # p-value = 0.3401
# Levene's test for homogeneity of variances   
leveneTest(size$beta[size$group %in% c("large", "small")& size$Site.x=="Center"] ~ size$group[size$group %in% c("large", "small")& size$Site.x=="Center"])     # 0.7383
# Step 2: Conduct the Parametric Test (t-test)
# Assuming assumptions are met
# Perform independent t-test
t.test(size$beta[size$group %in% c("large", "small")& size$Site.x=="Center"] ~ size$group[size$group %in% c("large", "small")& size$Site.x=="Center"])

shapiro.test(size$beta[size$group == "large" & size$Site.x=="Edge"])   # p-value = 0.3663
shapiro.test(size$beta[size$group == "small" & size$Site.x=="Edge"])   # p-value = 0.6378
# Levene's test for homogeneity of variances   
leveneTest(size$beta[size$group %in% c("large", "small")& size$Site.x=="Edge"] ~ size$group[size$group %in% c("large", "small")& size$Site.x=="Edge"])     # 0.1295
# Step 2: Conduct the Parametric Test (t-test)
# Assuming assumptions are met
# Perform independent t-test
t.test(size$beta[size$group %in% c("large", "small")& size$Site.x=="Edge"] ~ size$group[size$group %in% c("large", "small")& size$Site.x=="Edge"])

##########################
#Join together tables of each part...?
dist_beta <- dist_beta_0[,1:7]
i<-c(3, 4, 5, 6, 7)
dist_beta[ , i] <- apply(dist_beta[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))

#So it's dist_beta column two in these following items that are the problem...                          
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
acari_beta<-acari_beta[,c(7, 4, 5)]     # logdist, beta, Site.x
acari_beta$metric <- "3% OTU"
dist_beta<-dist_beta[,c(2, 3, 4)]       # logdist, beta, site
dist_beta$metric <- "zOTU"
names(dist_beta)<-c("log_dist", "beta", "Site.x", "metric")

beta <- rbind(dist_beta, acari_beta) 
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
                         
#######################################################
# Plot these
                         
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
  theme(       legend.position = "none", 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5),
        plot.title=element_text(size=70))
                         
b<- ggplot() + 
  geom_boxplot(data=beta[beta$metric=="3% OTU" & !is.na(beta$Site.x),],aes(x=Site.x, y=beta, fill=Site.x), color="black", size=1)+
  #facet_wrap(~metric, scales="free") +
  scale_fill_manual(values=SiteColors) +
  labs(title="B.", y="3% OTU beta diversity", x="") +
  KipukaTheme +
  guides(fill="none")+                       
  theme(strip.text = element_text(size = 70), 
        axis.text.x = element_text(angle=45, size=70, hjust=1, vjust=1), 
        axis.text.y = element_text(angle=45, size=70, margin=margin(0,-10,0,0)), 
        axis.title.y = element_text(size = 70, margin=margin(0,-25,0,0)),
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=70), 
        plot.title=element_text(size=70), 
        legend.text=element_text(size=60), 
        legend.title = element_blank(),
       legend.position = "none", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))   
                         
c <- ggplot() + 
  geom_point(data=nmds[nmds$Site!=c("C-E", "C-F"),],aes(x=MDS1zOTU,y=MDS2zOTU,colour=Site, size=pointsize, shape=Site), alpha=0.70, stroke=3) + 
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge", "Lava", "Kona", "Stainback")) +
  scale_shape_manual("Site", values=c("Center" = 16, "Edge" = 16, "Lava"=3, "Kona"=2, "Stainback"=2)) +
  scale_size_continuous("Kipuka area ("~m^2~")", range=c(2,32), breaks=seq(2,32,5), labels=round((10*seq(2,32,5))^2,100)) +
  labs(title="C.", x="NMDS1", y="NMDS2") +
  #coord_equal() +
  scale_x_continuous(breaks=seq(-2,1.5,0.5)) +
  guides(colour = guide_legend(override.aes = list(size=4))) + 
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
        plot.title=element_text(size=70))                     
                     
d<- ggplot() + 
  geom_boxplot(data=beta[beta$metric=="zOTU" & !is.na(beta$Site.x),],aes(x=Site.x, y=beta, fill=Site.x), color="black", size=1)+
  #facet_wrap(~metric, scales="free") +
  scale_fill_manual(values=SiteColors) +
  labs(title="D.", y="zOTU beta diversity", x="") +
  KipukaTheme +
  guides(fill="none")+                       
  theme(strip.text = element_text(size = 70), 
        axis.text.x = element_text(angle=45, size=70, hjust=1, vjust=1), 
        axis.text.y = element_text(angle=45, size=70, margin=margin(0,-10,0,0)), 
        axis.title.y = element_text(size = 70, margin=margin(0,-25,0,0)),
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=70), 
        plot.title=element_text(size=70), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top", 
        plot.margin = margin(0.2,1,0,1.35, "cm"))  
                         
jpeg("../Figures/NMDS-turnovers.jpg", width=2000, height=2000) 
plot_grid(a,b,c,d, nrow=2, ncol=2, rel_widths=c(1, 0.8), rel_heights=c(1,1))                         
dev.off()                                 

#################################################################################
# First, do a PERMANOVA to look for overall community comp differences between sites

# FIRST, 3% RADIUS OTU... ############################################                   
# curate the distance matrix                         
rownames(otu) <- gsub("X", "", rownames(otu))# Remove the character "X" from column and row names
colnames(otu) <- gsub("X", "", colnames(otu))
otu <- otu[complete.cases(otu), , drop = FALSE] # remove the part of otu below the distance matrix
rownames(otu) <- otu[[1]] # now make the first colum the row names. 
otu <- otu[, -1]                                                  
# Order the richness$Site data in the order of the names(otu),
# dropping out those smallest kipuka that have been removed from this analysis
Site <- richness[richness$ï..ID %in% names(otu), , drop = FALSE]
Site <- as.factor(Site$Site[order(match(Site$ï..ID, names(otu)))])                         
adonis2(otu ~ Site, permutations = 999)
                       
# Now, the pairwise test:
# First, I need a 3% radius OTU "OTU table"
vegan_otu <- function(otu_tab) {
  OTU <- otu_table(otu_tab, taxa_are_rows=F)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)}
  return(as(OTU, "matrix"))}
nmds<-read.csv("nmds3otu.csv")                    
my_otu_tab<-nmds[order(match(nmds$ID, names(otu))),]    
my_otu_tab<-as.data.frame(my_otu_tab[24:ncol(my_otu_tab)])                      
pairwise.adonis(vegan_otu(my_otu_tab), Site, p.adjust.m="fdr")
                         
# NEXT, ZOTU TESTS... #############################################
zOTU<-read.csv("zotu.csv")
zOTU <- zOTU[complete.cases(zOTU), , drop = FALSE] # remove the part of otu below the distance matrix
rownames(zOTU) <- gsub("X", "", rownames(zOTU))# Remove the character "X" from column and row names
colnames(zOTU) <- gsub("X", "", colnames(zOTU))
rownames(zOTU) <- zOTU[[1]] # now make the first colum the row names. 
zOTU <- zOTU[, -1]                                                 
# Order the richness$Site data in the order of the names(otu),
# dropping out those smallest kipuka that have been removed from this analysis
Site <- richness[richness$ï..ID %in% names(zOTU), ]
Site <- Site$Site[order(match(Site$ï..ID, names(zOTU)))]                        
adonis2(zOTU ~ Site, permutations = 999)  
                         
# Now the pairwise test   
my_otu_tab<-richness[order(match(richness$ï..ID, names(zOTU))),]    
my_otu_tab<-my_otu_tab[31:ncol(my_otu_tab)]                      
pairwise.adonis(vegan_otu(my_otu_tab), Site, p.adjust.m="bonferroni")
        
#######################################################
# Test difference between area types in zOTU and 3% radius OTU turnover

for (i in 1:length(unique(beta$metric))){  
        level<-unique(beta$metric)[i]
        print(level)
        test<-beta[beta$metric==level,]
        # First, test assumptions:  Assumes normal distribution of data and equal variances between groups.               
        # Check normality of residuals
        residuals <- lm(beta ~ Site.x, data = test)$residuals
        par(mar = c(1, 1, 1, 1))
        qqPlot(residuals, main = "Normal Q-Q Plot of Residuals")
        # Check homogeneity of variances
        p1<-leveneTest(beta ~ Site.x, data = test)[1,3]
        print(paste0("p-value for Levene's is ",p1))
        # Shapiro-Wilk test for normality (optional)
        p2<-shapiro.test(residuals)$p.value
        print(paste0("p-value for Shapiro-Wilk's is ",p2))
        if (p1 < 0.05 || p2 < 0.05){   # if these tests are significant, conduct a Kruskal-wallis test
              kruskal_result <- kruskal.test(beta ~ Site.x, data = test)
              cat("Kruskal-Wallis test results for", level, "\n")
              print(kruskal_result)
              pairwise_result <- pairwise.wilcox.test(test$beta, test$Site.x, p.adj = "bonferroni")
              cat("Pairwise Wilcoxon test results for", level, "\n")
              print(pairwise_result)
        }
        else { # Conduct the ANOVA                   
                anova_result <- aov(beta ~ Site.x, data = test)
                print(paste0("ANOVA test results for ", level))
                print(summary(anova_result)) 
              # Post hoc Tukey's HSD test
              tukey_result <- TukeyHSD(anova_result)
              cat("Tukey's HSD test results for", level, "\n")
              print(tukey_result)
        }
}                                                
                     
#######################################################
# Test difference between area size and zOTU and 3% radius OTU turnover

# First, test assumptions:  
# Fit linear regression model
gam_model <- gam(dist ~ s(log10(Area)), data = CE[CE$metric == "3% OTU", ])
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
ggplot(data = CE[CE$metric == "3% OTU", ], aes(x = Area, y = dist)) + 
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
