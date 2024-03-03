#Project: Kipukas               #
#Script: Size_Richness.R     #
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
library(car)
library(lmtest)

richness <- read.csv("merged_by_site_2.csv")
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainback"="#999933")
#Establish some themes up top to apply to all
KipukaTheme <- theme(axis.title=element_text(size=50), 
        axis.text = element_text(size=25, angle=50), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=40), 
        legend.key.height = unit(1, "cm"), 
        #legend.key.width = unit(1.5,"cm"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.title = element_text(size=40), 
        text = element_text(family = "serif"), 
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.spacing.y = unit(0.1,"cm")) 
                   
richness_mod_2 <- melt(richness, idvars = c("SiteID", "Site"), measure = c("SR", "SROTU"))
richness_mod_2 <- richness_mod_2[order(richness_mod_2$value, decreasing = TRUE),]  

# New facet label names for supp variable
supp.labs <- c("zOTU richness", "3% OTU richness")
names(supp.labs) <- c("SR", "SROTU")                   

#Reorder facets
richness_mod_2$variable <- factor(richness_mod_2$variable, levels = rev(c("SR", "SROTU")))                     

#################################################################################
# TEST:  ANOVA to check whether 3%OTU and the zOTU is different for each "area type" ( lava, edge, center, Stainback, Kona)

for (i in 1:length(unique(richness_mod_2$variable))){  
        level<-unique(richness_mod_2$variable)[i]
        print(level)
        test<-richness_mod_2[richness_mod_2$variable==level,]
        # First, test assumptions:  Assumes normal distribution of data and equal variances between groups.               
        # Check normality of residuals
        residuals <- lm(value ~ Site, data = test)$residuals
        par(mar = c(1, 1, 1, 1))  
        qqPlot(residuals, main = "Normal Q-Q Plot of Residuals")
        # Check homogeneity of variances
        p1<-leveneTest(value ~ Site, data = test)[1,3]
        print(paste0("p-value for Levene's is ",p1))
        # Shapiro-Wilk test for normality (optional)
        p2<-shapiro.test(residuals)$p.value
        print(paste0("p-value for Shapiro-Wilk's is ",p2))
        if (p1 < 0.05 || p2 < 0.05){   # if these tests are significant, conduct a Kruskal-wallis test
              kruskal_result <- kruskal.test(value ~ Site, data = test)
              cat("Kruskal-Wallis test results for", level, "\n")
              print(kruskal_result)
              pairwise_result <- pairwise.wilcox.test(test$value, test$Site, p.adj = "bonferroni")
              cat("Pairwise Wilcoxon test results for", level, "\n")
              print(pairwise_result)
        }
        else { # Conduct the ANOVA                   
                anova_result <- aov(value ~ Site, data = test)
                print(paste0("ANOVA test results for ", level))
                print(summary(anova_result))     
              # Post hoc Tukey's HSD test
              tukey_result <- TukeyHSD(anova_result)
              cat("Tukey's HSD test results for", level, "\n")
              print(tukey_result)
        }
}

###############################################################################################
# Linear regression for size vs. richness

for (j in 1:length(unique(richness_mod_2$variable))) {
        level<-unique(richness_mod_2$variable)[j]
        print(level)        
        for (i in 1:length(c("Center", "Edge"))){  
                type<-c("Center", "Edge")[i]
                print(type)
                test<-richness_mod_2[richness_mod_2$Site==type & richness_mod_2$variable==level,]
                # First, test assumptions:  
                # Fit linear regression model
                glm_model <- glm(value ~ log10(Area), data = test, family = poisson)
                # pseudo R2
                null_model <- glm(value ~ 1, data = test, family = poisson)
                McFadden_R2 <- 1 - (glm_model$deviance / glm(null_model)$deviance)
                print(paste0("McFadden r2 is ",McFadden_R2))                
                # 1. Linearity Check (Use Residuals vs. Fitted plot)
                par(mar = c(1, 1, 1, 1))
                par(mar = c(1, 1, 1, 1))  
                plot(glm_model, which = 1)                
                # 3. Homoscedasticity Check (Use Residuals vs. Fitted plot)
                plot(glm_model, which = 3)
                # 3. normality of residuals-- Q-Q plot -- for Poisson models, migt not perfectly follow normal distribution
                qqnorm(resid(glm_model))
                qqline(resid(glm_model))
                # 4. Overdispersion -- if the ratio is much larger than 1, there might be overdispersion
                df_resid <- df.residual(glm_model)
                dev_over_df <- deviance(glm_model) / df_resid
                print(paste0("overdispersion ratio is ",dev_over_df))
                # 5. Influence and outliers -- cook's distance                
                #infl <- influence.measures(glm_model)
                #plot(infl, which = "cook")
                # 6. goodness-of-fit tests, such as the Pearson or deviance goodness-of-fit tests. A high p-value suggests good fit.
                p<-pchisq(deviance(glm_model), df = df_resid, lower.tail = FALSE)
                print(paste0("goodness of fit p-val is ",p))

    #            if (DW < 0.05 || ST < 0.05){   # if these tests are significant, conduct a Kruskal-wallis test
     #                 kruskal_result <- kruskal.test(value ~ Site, data = test)
      #                cat("Kruskal-Wallis test results for", level, type, "\n")
       #               print(kruskal_result)
        #        }
         #       else { # Conduct the ANOVA                   
                        print(paste0("linear regression for ", level, type))
                        print(summary(glm_model))                  
             #   }
        }
}



# Print the summary of the linear regression model
summary(lm_model)


CenterzOTU <- lm(richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SR"]~richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SR"])
Center3otu<-lm(richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SROTU"]~richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SROTU"])
EdgezOTU<-lm(richness_mod_2$value[richness_mod_2$variable=="SR" & richness_mod_2$Site=="Edge"]~richness_mod_2$value[richness_mod_2$variable=="SR" & richness_mod_2$Site=="Edge"])
Edge3otu<-lm(richness_mod_2$value[richness_mod_2$variable=="SROTU" & richness_mod_2$Site=="Edge"]~richness_mod_2$value[richness_mod_2$variable=="SROTU" & richness_mod_2$Site=="Edge"])

m <- lm(y ~ x)
r2 = format(summary(m)$r.squared, digits = 3)))       

################################################################################################
# Plot it 

a <- ggplot() + 
  geom_boxplot(data=richness_mod_2,aes(x=reorder(Site, value), y=value, fill=Site), color="black", size=1)+
  facet_wrap(~variable, scales="free", 
             labeller = labeller(variable = supp.labs)) +
  scale_fill_manual(values=SiteColors) +
  labs(title="A.", x="") +
  KipukaTheme +
  theme(strip.text = element_text(size = 35), 
        axis.text.y = element_text(angle=45, size=45), 
        axis.text.x = element_text(angle=45, size=45, vjust=0.6), 
        axis.title.y = element_blank(), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1), 
      plot.margin = unit(c(0, 0, 0, 0), "cm"), 
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=50), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=45, hjust=0.4), 
        legend.title = element_blank(),
        legend.key.width = unit(0.6,"cm"), 
       legend.position = "top")              
                  

b<-ggplot() + 
  geom_smooth(method='lm', data=richness_mod_2[richness_mod_2$Site=="Center" | richness_mod_2$Site=="Edge",], aes(x=Area, y=value, colour=Site, fill=Site, linetype=variable), size=1, alpha=0.20)+  
  geom_point(data=richness_mod_2[richness_mod_2$Site=="Center" | richness_mod_2$Site=="Edge",],aes(x=Area, y=value, colour=Site, shape=variable), alpha=0.70, size=6, stroke = 3) + 
  scale_shape_manual("Site", values=c("SR" = 0, "SROTU"=15), labels=c("SR"="zOTU","SROTU"="3% OTU")) +
  scale_colour_manual(values=SiteColors) +
  scale_fill_manual(values=SiteColors)+ 
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +                 
  facet_wrap(~Site)+                 
  labs(title="B.", x="Kipuka area ("~m^2~")", y="OTU richness") +
  KipukaTheme +
  guides(color="none", fill="none", linetype="none") + 
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title.y=element_text(size=50, vjust=-0.5, margin=margin(r=10)), 
       axis.title.x=element_text(size=50, vjust=2, margin=margin(t=10)), 
        axis.text.y = element_text(size=45), 
        axis.text.x = element_text(size=45, vjust=1), 
        plot.title=element_text(size=50), 
        #legend.key.width = unit(7,"cm"), 
        legend.text=element_text(size=45, hjust=0.5), 
        legend.title = element_blank(),
       legend.position = "top")

jpeg("Figures/Figure3.jpg", width=2000, height=1000)
plot_grid(a, b, ncol = 2, rel_widths = c(1, 2))
dev.off()                     
 
                     

                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
##########################################
#old plot format                     
b <- ggplot() + 
  geom_smooth(method='lm', data=richness[richness$Site=="Center" | richness$Site=="Edge",], aes(x=Arealog, y=SR, colour=Site, fill=Site), size=1, alpha=0.20)+ #linetype=variable, 
  geom_point(data=richness[richness$Site=="Center" | richness$Site=="Edge",],aes(x=Arealog, y=SR, colour=Site), alpha=0.70, size=6, stroke = 3) + #, shape=variable
  geom_hline(yintercept=mean(richness$SR[richness$Site=="Kona"]), colour="#117733", lwd=3)+
  geom_hline(yintercept=mean(richness$SR[richness$Site=="Stainback"]), colour="#999933", lwd=3)+
  geom_hline(yintercept=mean(richness$SR[richness$Site=="Lava"]), colour="#888888", lwd=3, alpha=0.6)+
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge", "Kona", "Stainback")) +
  scale_fill_manual(values=SiteColors)+  
  facet_wrap(~Site)+                              
  #scale_linetype_discrete(values=c(2,5)) +
  labs(title="B.", x="Log area ("~m^2~")", y="zOTU richness") +
  KipukaTheme +
  guides(color="none", fill="none") +#shape="none", 
  theme(strip.text = element_text(size = 35), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=50), 
        axis.text.y = element_text(size=45), 
        axis.text.x = element_text(size=45, vjust=0.6),  
        plot.title=element_text(size=50), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top")    
                   
################################################################################
#haplotype richness within OTUs for Kipuka Centers and Kipuka edges

c <- ggplot() + 
  geom_smooth(method='lm', data=richness[richness$Site=="Center" | richness$Site=="Edge",], aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site, fill=Site), size=1, alpha=0.20)+ #linetype=variable, 
  geom_point(data=richness[richness$Site=="Center" | richness$Site=="Edge",],aes(x=Arealog, y=HaplotypeRichnessWithin, colour=Site), alpha=0.70, size=6, stroke = 3) + #, shape=variable
  geom_hline(yintercept=mean(richness$HaplotypeRichnessWithin[richness$Site=="Kona"]), colour="#117733", lwd=3)+
  geom_hline(yintercept=mean(richness$HaplotypeRichnessWithin[richness$Site=="Stainback"]), colour="#999933", lwd=3)+
  geom_hline(yintercept=mean(richness$HaplotypeRichnessWithin[richness$Site=="Lava"]), colour="#888888", lwd=3, alpha=0.6)+
  facet_wrap(~Site)+
  scale_colour_manual(values=SiteColors, limits = c("Center", "Edge", "Kona", "Stainback")) +
  scale_fill_manual(values=SiteColors, limits = c("Center", "Edge"))+                 
  labs(title="C.", x="Log area ("~m^2~")", y="Haplotype richness within OTUs") +
  KipukaTheme +
  guides(color="none", fill="none") +
  theme(strip.text = element_text(size = 30), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title=element_text(size=50), 
        axis.text = element_text(size=45), 
        plot.title=element_text(size=50), 
        legend.text=element_text(size=45), 
        legend.title = element_blank(),
       legend.position = "top")

###########################################################
#Plot these three together

#jpeg("Figures/Figure3.jpg", width=2000, height=1000)
#plot_grid(a, b, ncol = 2, rel_widths = c(1, 1))
#dev.off()
                   
