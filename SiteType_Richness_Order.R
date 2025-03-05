#Set up my working environment
setwd("H:/My Drive/Kipuka/Data")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
library(scales)
library(car)
library(lmtest)
library(dplyr)
library(vegan)
library(lme4)
library(MASS)
library(pscl)

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainback"="#999933")

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


richness <- read.csv("merged_by_site_2.csv")
OTUtoKeep<-as.data.frame(t(richness[c(12,17:nrow(richness)), 33:ncol(richness)]))
names(OTUtoKeep)[1:3]<- c("Order","OTU", "zOTU")
names(OTUtoKeep)[4:ncol(OTUtoKeep)] <- richness[19:nrow(richness),1]
OTUtoKeep <- OTUtoKeep[OTUtoKeep$Order %in% c("Araneae", "Lepidoptera", "Coleoptera", "Diptera", "Psocoptera", "Hemiptera"),]
# Exclude blank (NA or empty string) values and create a table for the 17th column
OTUtoKeep <- OTUtoKeep[grepl("OTU", OTUtoKeep[[2]]), ]

# Identify the values in column 1 that occur more than once
values_to_keep <- names(which(table(OTUtoKeep[[2]]) > 1))
# Subset the dataframe to keep only rows where column 1 matches those values
OTUtoKeep_filtered <- OTUtoKeep[OTUtoKeep[[2]] %in% values_to_keep, ]
OTUtoKeep_filtered <- OTUtoKeep_filtered %>%
  mutate(across(4:ncol(OTUtoKeep_filtered), as.numeric))
summary_data <- OTUtoKeep_filtered %>%
  # Ensure columns 3:ncol(OTUtoKeep_filtered) are numeric
  mutate(across(4:ncol(OTUtoKeep_filtered), as.numeric)) %>%
  # Transform values greater than 0 to 1
  mutate(across(4:ncol(OTUtoKeep_filtered), ~ ifelse(. > 0, 1, 0)))
names(summary_data)[2]<- c("OTU")

# OTU counts per site per taxon... 
siteOTU <- as.data.frame(OTUtoKeep %>%
  mutate(across(4:ncol(OTUtoKeep), as.numeric)) %>%
  mutate(across(4:ncol(OTUtoKeep), ~ ifelse(. > 0, 1, 0))) %>%
  group_by(OTU, Order) %>%  # Group by the first column
  summarise(across(2:(ncol(OTUtoKeep)-2), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>%  # Take max in each group 
  group_by(Order) %>%                       
  summarise(across(2:(ncol(OTUtoKeep)-2), sum, na.rm = TRUE)))  # Sum the max values across groups 
row.names(siteOTU) <- paste0(siteOTU$Order)
siteOTU <- siteOTU[-1]

# zOTU counts per site per taxo... 
sitezOTU <- as.data.frame(summary_data %>%
  mutate(across(4:ncol(summary_data), as.numeric)) %>%
  group_by(Order) %>%  
  summarise(across(3:54, sum, na.rm = TRUE)))
# OTU counts for zOTU weighting...  
zOTUotu <- as.data.frame(summary_data %>%
  mutate(across(4:ncol(summary_data), as.numeric)) %>%
  mutate(across(4:ncol(summary_data), ~ ifelse(. > 0, 1, 0))) %>%
  group_by(OTU, Order) %>%  # Group by the first column
  summarise(across(2:(ncol(summary_data)-2), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>%  # Take max in each group 
  group_by(Order) %>%                       
  summarise(across(2:(ncol(summary_data)-2), sum, na.rm = TRUE)))  # Sum the max values across groups 
weighted_zOTU <- sitezOTU[2:ncol(sitezOTU)] / zOTUotu[2:ncol(zOTUotu)]
row.names(weighted_zOTU) <- paste0(sitezOTU$Order)

richness_mod_2 <- as.data.frame(cbind(richness$X[19:nrow(richness)], richness$X.8[19:nrow(richness)], richness$X.9[19:nrow(richness)], t(siteOTU), t(weighted_zOTU)))
names(richness_mod_2)[1:3] <- c("my_ID", "Site", "Area")
richness_mod_2$Area<-as.numeric(gsub(",","",as.character(richness_mod_2$Area)))
richness_mod_2$Site<-gsub("Stainbeck","Stainback",as.character(richness_mod_2$Site))
richness_mod_2[, 4:ncol(richness_mod_2)] <- lapply(richness_mod_2[, 4:ncol(richness_mod_2)], as.numeric)

##########################################
# 3 % OTU
richness_mod_2 <- as.data.frame(cbind(richness$X[19:nrow(richness)], richness$X.8[19:nrow(richness)], richness$X.9[19:nrow(richness)], t(siteOTU), t(weighted_zOTU)))
names(richness_mod_2)[1:3] <- c("my_ID", "Site", "Area")
richness_mod_2$Area<-as.numeric(gsub(",","",as.character(richness_mod_2$Area)))
richness_mod_2$Site<-gsub("Stainbeck","Stainback",as.character(richness_mod_2$Site))
richness_mod_2[, 4:ncol(richness_mod_2)] <- lapply(richness_mod_2[, 4:ncol(richness_mod_2)], as.numeric)

#Let's just try all the orders! 
richness_mod_0 <- melt(richness_mod_2, idvars = c("SiteID", "Area"), measure = c("Araneae", "Psocoptera", "Hemiptera", "Lepidoptera", "Coleoptera", "Diptera"))
richness_mod_0$Area <- round(richness_mod_0$Area, 0)
richness_mod_0$Area[richness_mod_0$Site=="Lava" & is.na(richness_mod_0$Area)] <- "Lava"
richness_mod_0$Area[richness_mod_0$Site=="Kona" & is.na(richness_mod_0$Area)] <- "Kona"
richness_mod_0$Area[richness_mod_0$Site=="Stainbeck" & is.na(richness_mod_0$Area)] <- "Stainbeck"
richness_mod_0$Site <- factor(richness_mod_0$Site, levels=c("Lava", "Edge", "Center", "Stainback", "Kona"))

jpeg("../Figures/Order_Richness_OTU.jpg", width=3000, height=2000)
ggplot() + 
  geom_boxplot(data=richness_mod_0,aes(x=Site, y=value, fill=Site), color="black", size=1)+
  scale_fill_manual(values=SiteColors) +
  facet_wrap(~variable, nrow=2, ncol=3,scales="free_y") +
  guides(fill=guide_legend(nrow=2)) +
  labs(title="", x="", y="3% radius OTU richness") +
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
       legend.position = "none", 
        plot.margin = margin(1,1,.01,1, "cm"))
dev.off()

#################################################################################
# TEST:  ANOVA to check, for each taxon, whether zOTU richness is different for each "area type" ( lava, edge, center, Stainback, Kona)

for (i in 1:length(unique(richness_mod_0$variable))){  
        taxon<-unique(richness_mod_0$variable)[i]
        print(taxon)
        test<-richness_mod_0[richness_mod_0$variable==taxon,]
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
              cat("Kruskal-Wallis test results for", taxon, "\n")
              print(kruskal_result)
              pairwise_result <- pairwise.wilcox.test(test$value, test$Site, p.adj = "bonferroni")
              cat("Pairwise Wilcoxon test results for", taxon, "\n")
              print(pairwise_result)
        }
        else { # Conduct the ANOVA                   
                anova_result <- aov(value ~ Site, data = test)
                print(paste0("ANOVA test results for ", taxon))
                print(summary(anova_result))    
              # Post hoc Tukey's HSD test
              tukey_result <- TukeyHSD(anova_result)
              cat("Tukey's HSD test results for", taxon, "\n")
              print(tukey_result)
        }
}

################################################################################################

# REGRESS RICHNESS BY SIZE

# I need to obtain a dataframe with sites as rows, and id, type, area, 3% OTU, and zOTU as columns
# It would be helpful to iterate taxon-by-taxon

# Create a data version with weighted zOTU and 3% OTU
for (taxon in 1:nrow(weighted_zOTU)){
        richness_mod <- as.data.frame(cbind(richness$X[19:nrow(richness)], richness$X.8[19:nrow(richness)], richness$X.9[19:nrow(richness)], t(siteOTU[taxon,]), t(weighted_zOTU[taxon,])))
        names(richness_mod)[1:5] <- c("my_ID", "Site", "Area", "SROTU", "zOTU")
        richness_mod$Area<-as.numeric(gsub(",","",as.character(richness_mod$Area)))
        richness_mod[, 4:ncol(richness_mod)] <- lapply(richness_mod[, 4:ncol(richness_mod)], as.numeric)
        richness_mod <- richness_mod[richness_mod$Site %in% c("Center", "Edge"),]
        richness_mod$Area <- as.numeric(richness_mod$Area)
        richness_mod <- melt(richness_mod, idvars = c("Site", "Area"), measure = c("SROTU", "zOTU"))
        TAX<-row.names(weighted_zOTU)[taxon]
        print(TAX)
        for (i in 1:length(c("Center", "Edge"))){  
                type<-c("Center", "Edge")[i]
                print(type)
                test<-richness_mod[richness_mod$Site==type & richness_mod$variable=="SROTU",]
                # Fit linear regression model      
                model <- glm(value ~ log10(Area), data = test, family = poisson)
                null <- glm(value ~ 1, data = test, family = poisson)
                # Overdispersion -- if the ratio is much larger than 1, there might be overdispersion
                dispersion_ratio <- sum(residuals(glm_model, type = "pearson")^2) / df.residual(glm_model)
                if (dispersion_ratio>1.5) {
                        print(paste0("was over-dispersed... used negative binomial model for 3% OTU", TAX, type))
                        model <- glm.nb(value ~ log10(Area), data = test)
                        null <- glm.nb(value ~ 1, data = test)  # Null model (intercept only)                
                        # 1. Linearity Check (Use Residuals vs. Fitted plot)
                        par(mar = c(1, 1, 1, 1))
                        plot(model, which = 1)                
                        # 3. Homoscedasticity Check (Use Residuals vs. Fitted plot)
                        plot(model, which = 3)
                        # 3. normality of residuals-- Q-Q plot -- for Poisson models, migt not perfectly follow normal distribution
                        qqnorm(resid(model))
                        qqline(resid(model))
                        print(summary(model))
                        # Calculate adjusted McFadden's r2
                        adjusted_r2_nb <- 1 - (logLik(model) / logLik(null))
                        print(paste0("mcfadden's r2 is ", adjusted_r2_nb))
                        }
                else if (dispersion_ratio<0.5) {
                        print(paste0("was under-dispersed... used gam for 3% OTU", TAX, type))
                        model <- gam(value ~ s(log10(Area)), family = poisson, data = test)
                        null <- gam(value ~ 1, family = poisson, data = test)  # Null model (intercept only)                
                        print(summary(model))
                        # 1. Linearity Check (Use Residuals vs. Fitted plot)
                        plot(model, residuals = TRUE, pch = 20, main = "Smoothed Relationship: log10(Area)")                
                        # 3. Homoscedasticity Check (Use Residuals vs. Fitted plot)
                        plot(fitted(model), residuals(gam_model), main = "Residuals vs Fitted", xlab = "Fitted Values", ylab = "Residuals", pch = 20)
                        abline(h = 0, col = "red", lty = 2)
                        # 3. normality of residuals-- Q-Q plot -- for Poisson models, migt not perfectly follow normal distribution
                        qqnorm(residuals(model))
                        qqline(residuals(model), col = "red") 
                        # Calculate adjusted McFadden's r2
                        logLik_gam <- logLik(model)[1]  # Extract log-likelihood
                        logLik_null <- logLik(null)[1]  # Extract log-likelihood of null model
                        k_gam <- model$df.residual  # Effective degrees of freedom
                        # Compute Adjusted McFadden's R²
                        adjusted_r2_gam <- 1 - (logLik(model) / logLik(null))
                        print(paste0("mcfadden's r2 is: ", round(adjusted_r2_gam, 4)))
                }
                else {               
                # 1. Linearity Check (Use Residuals vs. Fitted plot)
                # points should be randomly scattered above/below
                par(mar = c(1, 1, 1, 1))
                plot(model, which = 1)                
                # 2. Homoscedasticity Check (Use Residuals vs. Fitted plot)
                # no funnel or cone shape should be visible here
                plot(model, which = 3)
                # 3. normality of residuals-- Q-Q plot -- for Poisson models, migt not perfectly follow normal distribution
                # points should fall along the line
                qqnorm(resid(model))
                qqline(resid(model))
                # 5. Influence and outliers -- cook's distance                
                infl <- influence.measures(model)
                cooks_d <- cooks.distance(model)
                # Plot Cook's distances
                plot(cooks_d, type = "h", main = "Cook's Distances", ylab = "Cook's Distance", xlab = "Observation Index",col = "blue")
                # Add a reference line for a threshold (commonly 4/n)
                abline(h = 4 / nrow(model$model), col = "red", lty = 2)                
                which(cooks_d > (4 / nrow(model$model)))                
                print(paste0("linear regression for 3% OTU", TAX, type))
                print(summary(model)) 
                #McFadden_R2 <- 1 - (glm_model$deviance / glm(null_model)$deviance)
                #print(McFadden_adj_R2)
                #McFadden_adj_R2 <- 1 - (logLik(model) / logLik(null))
                #print(paste0("mcfadden's r2 is ",McFadden_adj_R2))
                print(pR2(model))
                }
                
                #Starting with zOTUs
                #print("zOTU")
                #test<-richness_mod[richness_mod$Site==type & richness_mod$variable=="zOTU",]
                # Fit linear regression model. Because zOTU is now continuous, we can't use poisson...       
                #lm_model <- lm(value ~ log10(Area), data = test)  # Equivalent to Gaussian GLM
                # 1. Linearity Check (Use Residuals vs. Fitted plot)
                # points should be randomly scattered above/below
                #par(mar = c(1, 1, 1, 1))
                #plot(lm_model, which = 1)                
                # 2. Homoscedasticity Check (Use Residuals vs. Fitted plot)
                # no funnel or cone shape should be visible here
                #plot(lm_model, which = 3)
                # 3. normality of residuals-- Q-Q plot -- for Poisson models, migt not perfectly follow normal distribution
                # points should fall along the line
                #qqnorm(resid(lm_model))
                #qqline(resid(lm_model))
                # 5. Influence and outliers -- cook's distance                
                #infl <- influence.measures(lm_model)
                #cooks_d <- cooks.distance(lm_model)
                # Plot Cook's distances
                #plot(cooks_d, type = "h", main = "Cook's Distances", ylab = "Cook's Distance", xlab = "Observation Index",col = "blue")
                # Add a reference line for a threshold (commonly 4/n)
                #abline(h = 4 / nrow(lm_model$model), col = "red", lty = 2)                
                #which(cooks_d > (4 / nrow(lm_model$model)))
                #print(summary(lm_model)) 
        }
}


                                                            



#Put in correct order
#richness_mod_0$Area <- factor(richness_mod_0$Area, levels=c("Lava", "3", "4", "5", "Stainbeck", "Kona"))

richness_mod_0$Area<-as.numeric(richness_mod_0$Area)
richness_mod_0$linetype<-"dashed"
richness_mod_0$linetype[richness_mod_0$variable %in% c("Coleoptera", "Diptera", "Lepidoptera", "Psocoptera") & richness_mod_0$Site=="Center"]<-"solid"
richness_mod_0$linetype[richness_mod_0$variable %in% c("Araneae", "Hemiptera") & richness_mod_0$Site=="Edge"]<-"solid"

jpeg("../Figures/Order_Richness_OTU_size.jpg", width=3000, height=2000)
ggplot() + 
  geom_point(data=richness_mod_0,aes(x=Area, y=value, colour=Site, fill=Site), size=4.5, shape=0, stroke=2)+
  geom_smooth(method='lm', data=richness_mod_0, aes(x=Area, y=value, colour=Site, fill=Site, alpha=linetype, linetype=linetype))+
  scale_colour_manual(values=SiteColors) +
  facet_wrap(~variable, nrow=2, scales="free") +
  scale_fill_manual(values = SiteColors) + 
  guides(fill=guide_legend(nrow=1), colour="none") +
  scale_linetype_manual(values = c("dashed"=2, "solid"=1)) +
  scale_alpha_manual(values = c("dashed"=0, "solid"=.4))+ 
  labs(title="", x="Kipuka area ("~m^2~")", y="3% radius OTU richness") +
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  + 
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
       axis.title.x=element_text(size=55, margin=margin(-15,0,0,0)), 
       axis.title.y=element_text(size=55, margin=margin(0,-15,0,0)), 
        axis.text.y = element_text(size=50, angle=45), 
        axis.text.x = element_text(size=50, angle=45, vjust=0.5, hjust=0.5), 
        plot.title=element_text(size=55), 
        legend.text=element_text(size=50), 
        legend.title = element_text(size=50),
       legend.position = "none", 
        plot.margin = margin(1,1,.01,1, "cm"))
dev.off()












##########################################
# zOTU

richness_mod_2 <- as.data.frame(cbind(richness$X[19:nrow(richness)], richness$X.8[19:nrow(richness)], richness$X.9[19:nrow(richness)], t(weighted_zOTU)))
names(richness_mod_2)[1:3] <- c("my_ID", "Site", "Area")
richness_mod_2$Area<-as.numeric(gsub(",","",as.character(richness_mod_2$Area)))
richness_mod_2$Site<-gsub("Stainbeck","Stainback",as.character(richness_mod_2$Site))
richness_mod_2[, 4:ncol(richness_mod_2)] <- lapply(richness_mod_2[, 4:ncol(richness_mod_2)], as.numeric)

#Let's just try all the orders! 
richness_mod_0 <- melt(richness_mod_2, idvars = c("SiteID", "Area"), measure = c("Araneae", "Psocoptera", "Hemiptera", "Lepidoptera", "Coleoptera", "Diptera"))
richness_mod_0$Area <- round(richness_mod_0$Area, 0)
richness_mod_0$Area[richness_mod_0$Site=="Lava" & is.na(richness_mod_0$Area)] <- "Lava"
richness_mod_0$Area[richness_mod_0$Site=="Kona" & is.na(richness_mod_0$Area)] <- "Kona"
richness_mod_0$Area[richness_mod_0$Site=="Stainbeck" & is.na(richness_mod_0$Area)] <- "Stainbeck"
richness_mod_0$Site <- factor(richness_mod_0$Site, levels=c("Lava", "Edge", "Center", "Stainback", "Kona"))

jpeg("../Figures/Order_Richness_zOTU.jpg", width=3000, height=2000)
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
       legend.position = "none", 
        plot.margin = margin(1,1,.01,1, "cm"))
dev.off()
