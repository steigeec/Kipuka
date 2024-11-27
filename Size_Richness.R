
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


# REMOVE LEGEND TO THE BOTTOM OF A and B
# Make the lines IN b THICKER 


richness <- read.csv("merged_by_site_2.csv")
OTUtoKeep<-as.data.frame(t(richness[17:nrow(richness), 33:ncol(richness)]))
names(OTUtoKeep)[1]<- c("OTU")

# Exclude blank (NA or empty string) values and create a table for the 17th column
OTUtoKeep <- OTUtoKeep[grepl("OTU", OTUtoKeep[[1]]), ]
# Identify the values in column 1 that occur more than once
values_to_keep <- names(which(table(OTUtoKeep[[1]]) > 1))

# Subset the dataframe to keep only rows where column 1 matches those values
OTUtoKeep_filtered <- OTUtoKeep[OTUtoKeep[[1]] %in% values_to_keep, ]

OTUtoKeep_filtered <- OTUtoKeep_filtered %>%
  mutate(across(3:ncol(OTUtoKeep_filtered), as.numeric))

summary_data <- OTUtoKeep_filtered %>%
  # Ensure columns 3:ncol(OTUtoKeep_filtered) are numeric
  mutate(across(3:ncol(OTUtoKeep_filtered), as.numeric)) %>%
  # Transform values greater than 0 to 1
  mutate(across(3:ncol(OTUtoKeep_filtered), ~ ifelse(. > 0, 1, 0)))
names(summary_data)[1]<- c("OTU")

# OTU counts per site...
siteOTU <- as.list(OTUtoKeep %>%
  mutate(across(3:ncol(OTUtoKeep), as.numeric)) %>%
  mutate(across(3:ncol(OTUtoKeep), ~ ifelse(. > 0, 1, 0))) %>%                
  group_by(OTU) %>%  # Group by the first column
  summarise(across(2:(ncol(summary_data)-1), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>%  # Take max in each group 
  summarise(across(2:(ncol(summary_data)-1), sum, na.rm = TRUE)))  # Sum the max values across groups 

# zOTU counts per site... 
sitezOTU <- as.list(summary_data %>%
  summarise(across(3:ncol(summary_data), sum, na.rm = TRUE)))
# OTU counts for zOTU weighting...  
zOTUotu <- as.list(summary_data %>%
  group_by(OTU) %>%  # Group by the first column
  summarise(across(2:(ncol(summary_data)-1), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>%  # Take max in each group 
  summarise(across(2:(ncol(summary_data)-1), sum, na.rm = TRUE)))  # Sum the max values across groups 

richness_mod_2 <- as.data.frame(cbind(richness$X[19:nrow(richness)], richness$X.8[19:nrow(richness)], richness$X.9[19:nrow(richness)], siteOTU, sitezOTU, zOTUotu))
names(richness_mod_2) <- c("my_ID", "Site", "Area", "SROTU", "unweighted_zOTU", "OTU")
richness_mod_2$Area<-as.numeric(gsub(",","",as.character(richness_mod_2$Area)))
richness_mod_2$Site<-gsub("Stainbeck","Stainback",as.character(richness_mod_2$Site))
richness_mod_2[, 3:6] <- lapply(richness_mod_2[, 3:6], as.numeric)
#Weight zoTU by OTU richness per site
richness_mod_2$zOTU<-richness_mod_2$unweighted_zOTU/richness_mod_2$OTU

richness_mod_2 <- melt(richness_mod_2, idvars = c("my_ID", "Site","Area", "OTU", "unweighted_zOTU"), measure = c("zOTU", "SROTU"))
richness_mod_2 <- richness_mod_2[order(richness_mod_2$value, decreasing = TRUE),]  

# New facet label names for supp variable
supp.labs <- c("zOTU richness", "3% OTU richness")
names(supp.labs) <- c("zOTU", "SROTU")                   

#Reorder facets
richness_mod_2$variable <- factor(richness_mod_2$variable, levels = rev(c("zOTU", "SROTU")))   
richness_mod_2$Site <- factor(richness_mod_2$Site, levels = unique(richness_mod_2$Site))
richness_mod_2$my_ID <- factor(richness_mod_2$my_ID, levels = unique(richness_mod_2$my_ID))
richness_mod_2 <- richness_mod_2 %>%
  mutate(alpha_value = ifelse(Site == "Center" & variable == "SROTU", 0.2, 0))
richness_mod_2 <- richness_mod_2 %>%
  mutate(group = case_when(
    Site %in% c("Kona", "Stainback") ~ "forest",
    Site %in% c("Center", "Edge") ~ "kipuka",
    Site %in% c("Lava") ~ "lava"  ))

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




# Fit a linear mixed model (example)
model <- lmer(value ~ group + (1 | Site), data = richness_mod_2[richness_mod_2$variable=="SROTU",])

# Check summary
summary(model)



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
                # points should be randomly scattered above/below
                par(mar = c(1, 1, 1, 1))
                par(mar = c(1, 1, 1, 1))  
                plot(glm_model, which = 1)                
                # 2. Homoscedasticity Check (Use Residuals vs. Fitted plot)
                # no funnel or cone shape should be visible here
                plot(glm_model, which = 3)
                # 3. normality of residuals-- Q-Q plot -- for Poisson models, migt not perfectly follow normal distribution
                # points should fall along the line
                qqnorm(resid(glm_model))
                qqline(resid(glm_model))
                # 4. Overdispersion -- if the ratio is much larger than 1, there might be overdispersion
                df_resid <- df.residual(glm_model)
                dev_over_df <- deviance(glm_model) / df_resid
                print(paste0("overdispersion ratio is ",dev_over_df))
                # 5. Influence and outliers -- cook's distance                
                infl <- influence.measures(glm_model)
                cooks_d <- cooks.distance(glm_model)
                # Plot Cook's distances
                plot(cooks_d, type = "h", 
                     main = "Cook's Distances", 
                     ylab = "Cook's Distance", 
                     xlab = "Observation Index",
                     col = "blue")
                # Add a reference line for a threshold (commonly 4/n)
                abline(h = 4 / nrow(glm_model$model), col = "red", lty = 2)                
                which(cooks_d > (4 / nrow(glm_model$model)))
                # 6. goodness-of-fit tests, such as the Pearson or deviance goodness-of-fit tests. A high p-value suggests good fit.
                p<-pchisq(deviance(glm_model), df = df_resid, lower.tail = FALSE)
                print(paste0("goodness of fit p-val is ",p))                  
                print(paste0("linear regression for ", level, type))
                print(summary(glm_model))                  
        }
}

# Print the summary of the linear regression model

CenterzOTU <- lm(richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="zOTU"]~richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="zOTU"])
Center3otu<-lm(richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SROTU"]~richness_mod_2$value[richness_mod_2$Site=="Center" & richness_mod_2$variable=="SROTU"])
EdgezOTU<-lm(richness_mod_2$value[richness_mod_2$variable=="zOTU" & richness_mod_2$Site=="Edge"]~richness_mod_2$value[richness_mod_2$variable=="zOTU" & richness_mod_2$Site=="Edge"])
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
       legend.position = "none")                           

b<-ggplot() + 
  geom_smooth(method='lm', 
              data = richness_mod_2[richness_mod_2$Site == "Center",], 
              aes(x = Area, y = value, colour = Site, fill = Site, linetype = variable, alpha = alpha_value), 
              size = 1) +  # Set a default alpha of 0.2 for smoothing line
  geom_point(data = richness_mod_2[richness_mod_2$Site == "Center",], 
             aes(x = Area, y = value, colour = Site, shape = variable), 
             size = 6, stroke = 3) +  # Use alpha_value for points only
  scale_shape_manual("Site", values = c("zOTU" = 0, "SROTU" = 15), labels = c("zOTU" = "zOTU", "SROTU" = "3% OTU")) +
  scale_colour_manual(values = SiteColors) +
  scale_fill_manual(values = SiteColors) + 
  scale_alpha_identity() +  # Use alpha_identity to apply the alpha value directly
  scale_x_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))  +                 
  facet_wrap(~variable, scales = "free", 
             labeller = labeller(variable = supp.labs)) +                 
  labs(title = "B.", x = "Kipuka area ("~m^2~")", y = "OTU richness") +
  KipukaTheme +
  guides(color = "none", fill = "none", linetype = "none") + 
  theme(strip.text = element_text(size = 45), 
        panel.grid.major = element_line(
          rgb(105, 105, 105, maxColorValue = 255),
          linetype = "dotted", 
          size = 1),   
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.minor = element_line(
          rgb(105, 105, 105, maxColorValue = 255),
          linetype = "dotted", 
          size = 0.5), 
        axis.title.y = element_text(size = 50, vjust = -0.5, margin = margin(r = 10)), 
        axis.title.x = element_text(size = 50, vjust = 2, margin = margin(t = 10)), 
        axis.text.y = element_text(size = 45), 
        axis.text.x = element_text(size = 45, vjust = 1), 
        plot.title = element_text(size = 50), 
        legend.text = element_text(size = 45, hjust = 0.5), 
        legend.title = element_blank(),
        legend.position = "none")
                 

jpeg("../Figures/Figure3.jpg", width=2000, height=1000)
plot_grid(a, b, ncol = 2, rel_widths = c(1, 1.1))
dev.off()                     
