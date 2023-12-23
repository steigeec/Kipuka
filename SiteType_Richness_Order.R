#Project: Kipukas                  #
#Script: SiteType_Richness_Order.R #
#Author: Emma Steigerwald          #
#Date:17 Mar 2022                  #
####################################

#Set up my working environment
setwd("G:/My Drive/Kipuka/Data")
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

#Let's just try all the orders! 
richness_mod_0 <- melt(richness, idvars = c("SiteID", "Arealog"), measure = c("Araneae", "Pscoptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Acari", "Coleoptera", "Diptera"))
richness_mod_0$Arealog <- round(richness_mod_0$Arealog, 0)
richness_mod_0$Arealog[richness_mod_0$Site=="Lava" & is.na(richness_mod_0$Arealog)] <- "Lava"
richness_mod_0$Arealog[richness_mod_0$Site=="Kona" & is.na(richness_mod_0$Arealog)] <- "Kona"
richness_mod_0$Arealog[richness_mod_0$Site=="Stainbeck" & is.na(richness_mod_0$Arealog)] <- "Stainbeck"

#Remove orders that aren't speciose enough
my_orders<-c("Araneae", "Coleoptera", "Diptera", "Hemiptera", "Lepidoptera", "Pscoptera")
richness_mod_0<-richness_mod_0[richness_mod_0$variable %in% my_orders,]

#Fix spelling error on sheet before proceeding
richness_mod_0$Site<-gsub("Stainbeck","Stainback",as.character(richness_mod_0$Site))

richness_mod_0$Site <- factor(richness_mod_0$Site, levels=c("Lava", "Edge", "Center", "Stainback", "Kona"))

#################################################################################
# TEST:  ANOVA to check, for each taxon, whether zOTU richness is different for each "area type" ( lava, edge, center, Stainback, Kona)

for (i in 1:length(unique(richness_mod_0$variable))){  
        taxon<-unique(richness_mod_0$variable)[i]
        print(taxon)
        test<-richness_mod_0[richness_mod_0$variable==taxon,]
        # First, test assumptions:  Assumes normal distribution of data and equal variances between groups.               
        # Check normality of residuals
        residuals <- lm(value ~ Site, data = test)$residuals
        qqPlot(residuals, main = "Normal Q-Q Plot of Residuals")
        # Check homogeneity of variances
        p1<-leveneTest(value ~ Site, data = test)[1,3]
        print(paste0("p-value for Levene's is ",p1))
        # Shapiro-Wilk test for normality (optional)
        p2<-shapiro.test(residuals)$p.value
        print(paste0("p-value for Shapiro-Wilk's is ",p2))
        if (p1 < 0.05 || p2 < 0.05){   # if these tests are significant, conduct a Kruskal-wallis test
              kruskal_result <- kruskal.test(value ~ Site, data = test)
              cat("Kruskal-Wallis test results for", taxon, "\n")
              print(kruskal_result)
        }
        else { # Conduct the ANOVA                   
                anova_result <- aov(value ~ Site, data = test)
                print(paste0("ANOVA test results for ", taxon))
                print(summary(anova_result))                  
        }
}

###################################################################################
# Now plot

jpeg("Figures/Order_Richness_1.jpg", width=3000, height=2000)
ggplot() + 
  geom_boxplot(data=richness_mod_0,aes(x=Site, y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  facet_wrap(~variable, nrow=2, ncol=3,scales="free_y") +
  guides(fill=guide_legend(nrow=2)) +
  labs(title="", x="", y="zOTU richness") +
  KipukaTheme +
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=55), 
        axis.text.y = element_text(size=50, angle=45), 
        axis.text.x = element_text(size=50, angle=45, vjust=0.5), 
        plot.title=element_text(size=55), 
        legend.text=element_text(size=50), 
        legend.title = element_text(size=50),
       legend.position = "top", 
        plot.margin = margin(1,1,.01,1, "cm"))
dev.off()
