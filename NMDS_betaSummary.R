
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

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainback"="#999933", "C-E"="black", "C-F"="white")
ExtendedSiteColors <- c("Edge_Center"="grey24", "Lava_Center"="grey24", "Stainback_Center"="grey24", "Kona_Center"="grey24", "Center_Edge"="grey24", 
                        "Lava_Edge"="grey24", "Stainback_Edge"="grey24", "Kona_Edge"="grey24", "Center"="#332288", "Edge"="#6699CC", 
                        "Lava"="#888888", "Stainback_Lava"="grey24", "Kona_Lava"="grey24", "Stainback"="#999933", "Kona_Stainback"="grey24", "Kona"="#117733")


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

# Switch out for jaccard... 
# Do I need to weight zoTU by OTU...?

richness <- read.csv("merged_by_site_2.csv")
names(richness) <- richness[18,]
richness <- richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))
richness$Site<-gsub("Stainbeck","Stainback",as.character(richness$Site))
#dist_beta_0 <- read.csv("Distance_v_beta.csv")
# beta diversity for 3 % OTU for within-area
otu <- read.csv("3OTU.csv")
# beta diversity for zOTU for within-area
zOTUbeta <- read.csv("zOTU_Bray.csv")
# geographic distances
geo_dist<-read.csv("geo_dist.csv")
nmds<-read.csv("nmds3otu.csv")
#Fix spelling error on sheet before proceeding
nmds$Site<-gsub("Stainbeck","Stainback",as.character(nmds$Site))
#BC<-read.csv("BC.csv")
# beta diversity for 3 % OTU for between-area
#BC3<-read.csv("BC3.csv")
# beta diversity for z OTU for between-area
#CE<-read.csv("C-E.csv")

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
# Produce a table of 3% OTU beta diversity values within and between habitat types
otu<-otu[complete.cases(otu),]
OTU3beta<-as.matrix(otu)
colnames(OTU3beta)<- gsub('[X]', '', colnames(OTU3beta))
rownames(OTU3beta)<-OTU3beta[,1]
OTU3beta<-OTU3beta[,-1]
OTU3beta<-data.frame(col=colnames(OTU3beta)[col(OTU3beta)], row=rownames(OTU3beta)[row(OTU3beta)], beta=c(OTU3beta))
#Add attributes of each site
#First we add whether it's center, edge, etc etc etc
OTU3beta<-merge(OTU3beta, richness[,c(1,9)], by.x="col", by.y="ID")
OTU3beta<-merge(OTU3beta, richness[,c(1,9)], by.x="row", by.y="ID")
#Now add the distances between these sites
OTU3beta$index<-paste(OTU3beta$row, OTU3beta$col, sep="_")
OTU3beta<-merge(OTU3beta, geo_dist[,3:5], by.x="index", by.y="index")
#remove the same-site pairs
OTU3beta<-OTU3beta[OTU3beta$row!=OTU3beta$col,]     
#Remove flipped pairs
OTU3beta <- OTU3beta%>% distinct(beta, .keep_all= TRUE)
OTU3beta$beta<-as.numeric(OTU3beta$beta)

# Now, for zOTU, weight zOTU per 3% OTU... 
zOTUbeta$index <- paste(zOTUbeta$Var1, zOTUbeta$Var2, sep="_")
colnames(zOTUbeta)[3] <- "zOTU"
zOTUbeta <- merge(zOTUbeta[,3:4], OTU3beta, by="index")

# Join on size data, because we will only use kipuka > 5000 m^2
zOTUbeta <- merge (zOTUbeta, richness[,c(1,10)], by.x="row", by.y="ID", all.x=T)
zOTUbeta <- merge (zOTUbeta, richness[,c(1,10)], by.x="col", by.y="ID", all.x=T)
zOTUbeta <- zOTUbeta[
  (is.na(zOTUbeta$Area.x) & is.na(zOTUbeta$Area.y)) | 
  (zOTUbeta$Area.x > 5000 & zOTUbeta$Area.y > 5000), 
]

# We need a category for coloring our box and whisker plots... 
zOTUbeta$Site<-paste(zOTUbeta$Site.x, zOTUbeta$Site.y, sep="_")
zOTUbeta$Site <- gsub("Center_Center", "Center", zOTUbeta$Site)
zOTUbeta$Site <- gsub("Edge_Edge", "Edge", zOTUbeta$Site)
zOTUbeta$Site <- gsub("Lava_Lava", "Lava", zOTUbeta$Site)
zOTUbeta$Site <- gsub("Stainback_Stainback", "Stainback", zOTUbeta$Site)
zOTUbeta$Site <- gsub("Kona_Kona", "Kona", zOTUbeta$Site)
zOTUbeta$Site <- factor(zOTUbeta$Site, levels=c("Lava", "Edge", "Center", "Stainback", "Kona", 
                                          "Edge_Center", "Lava_Center", "Stainback_Center", "Kona_Center", "Center_Edge", 
                                          "Lava_Edge", "Stainback_Edge", "Kona_Edge", "Stainback_Lava", "Kona_Lava", "Kona_Stainback"))  

#######################################################
# Plot these

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
  theme(axis.title.x = element_text(size=70), 
        axis.text.x = element_text(angle=45, size=70, hjust=1, vjust=1),
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
  geom_boxplot(data=zOTUbeta,aes(x=Site, y=beta, fill=Site), color="black", size=1)+
  #facet_wrap(~metric, scales="free") +
  scale_fill_manual(values=ExtendedSiteColors) +
  labs(title="B.", y="3% OTU beta diversity", x="") +
  KipukaTheme +
  guides(fill="none")+                       
  theme(axis.title.x=element_blank(), 
        strip.text = element_text(size = 70), 
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
    plot.title = element_text(size = 70, hjust=-.1), # Offset title to the left
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
        plot.title=element_text(size=70, hjust=-.1))                     
                     
d<- ggplot() + 
  geom_boxplot(data=zOTUbeta,aes(x=Site, y=zOTU, fill=Site), color="black", size=1)+
  scale_fill_manual(values=ExtendedSiteColors) +
  labs(title="C.", y="zOTU beta diversity", x="") +
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
        plot.title=element_text(size=70, hjust=-.1), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top", 
        #plot.margin = margin(0.2,1,0,.8, "cm")
       )  
                         
jpeg("../Figures/NMDS-turnovers.jpg", width=4500, height=1500) 
plot_grid(a,b,d, nrow=1, ncol=3, rel_widths=c(1.75, 2, 2))                         
dev.off()                                 

#######################################################
# Test for differences between

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

#################################################################################


















##########################
# With the table of within-area-type dissimilarity measures,  compare the dissimilarity within small edges to that within big edges,
# and the dissimilarity of small cores to that within big cores...

# Join the size data onto OTU3beta
size<-merge(OTU3beta, richness[c(1,10)], by.x="row", by.y="ï..ID")
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
