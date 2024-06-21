
#Set up my working environment
setwd("G:/My Drive/Kipuka")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
font_import()
library(vegan)
if (!requireNamespace("car", quietly = TRUE)) { install.packages("car")}
library(car)

richness <- read.csv("Data/merged_by_site_2.csv")
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))        

rep <- richness[1:22]
i<-seq(15, 21, 1)
rep[ , i] <- apply(rep[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))                         
rep$totalRichness <- rep$Hemiptera + rep$Lepidoptera + rep$Pscoptera +  rep$Araneae + rep$Coleoptera + rep$Diptera
                   
rep <- melt(rep, idvars = c("SiteID", "Site", "totalRichness"), measure.vars = c("Araneae", "Pscoptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Acari", "Coleoptera", "Diptera"))
rep$my_site <- paste(rep$Site, rep$SiteID)
rep$prop <- rep$value/rep$totalRichness
                   
#Fix spelling error on sheet before proceeding
rep$Site<-gsub("Stainbeck","Stainback",as.character(rep$Site))  

#I want ordered by my sites                   
rep$Site <- factor(rep$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava")))  
rep <- rename(rep, id = ï..ID) 
rep$Arealog<-round(as.numeric(rep$Arealog),1)
rep<-rep %>% dplyr::mutate(Arealog = tidyr::replace_na(Arealog, ""))                               
                   
#################################################################################
# TEST:  ANOVA to check, for each taxon, whether prop is different for each "area type" ( lava, edge, center, Stainback, Kona)

for (i in 1:length(unique(rep$variable))){  
        taxon<-unique(rep$variable)[i]
        print(taxon)
        test<-rep[rep$variable==taxon,]
        # First, test assumptions:  Assumes normal distribution of data and equal variances between groups.               
        # Check normality of residuals
        residuals <- lm(prop ~ Site, data = test)$residuals
        par(mar = c(1, 1, 1, 1))  
        qqPlot(residuals, main = "Normal Q-Q Plot of Residuals")
        # Check homogeneity of variances
        p1<-leveneTest(prop ~ Site, data = test)[1,3]
        print(paste0("p-value for Levene's is ",p1))
        # Shapiro-Wilk test for normality (optional)
        p2<-shapiro.test(residuals)$p.value
        print(paste0("p-value for Shapiro-Wilk's is ",p2))
        if (p1 < 0.05 || p2 < 0.05){   # if these tests are significant, conduct a Kruskal-wallis test
              kruskal_result <- kruskal.test(prop ~ Site, data = test)
              cat("Kruskal-Wallis test results for", taxon, "\n")
              print(kruskal_result)
              pairwise_result <- pairwise.wilcox.test(test$prop, test$Site, p.adj = "bonferroni")
              cat("Pairwise Wilcoxon test results for", taxon, "\n")
              print(pairwise_result)
        }
        else { # Conduct the ANOVA                   
                anova_result <- aov(prop ~ Site, data = test)
                print(paste0("ANOVA test results for ", taxon))
                print(summary(anova_result))       
              # Post hoc Tukey's HSD test
              tukey_result <- TukeyHSD(anova_result)
              cat("Tukey's HSD test results for", taxon, "\n")
              print(tukey_result)
        }
}
#################################################################################
# TEST:  ANOVA to check, for each taxon, whether value is different for each "area type" ( lava, edge, center, Stainback, Kona)

for (i in 1:length(unique(rep$variable))){  
        taxon<-unique(rep$variable)[i]
        print(taxon)
        test<-rep[rep$variable==taxon,]
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
              pairwise_result <- pairwise.wilcox.test(test$prop, test$Site, p.adj = "bonferroni")
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
                   
#################################################################################
#Proportional representation of orders by site                         
                   
#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainbeck"="#999933")

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
                   
a<- ggplot(data=rep, aes(x=reorder(my_site, Area), y=prop, width=1, fill=variable)) +
  geom_bar(stat="identity", color="black") + #, size=0.4, key_glyph = "polygon3"
  labs(title="A.") +
  xlab(expression("         [                 By increasing size ("~m^2~")                         ]                                                             "))+                 
  facet_grid(cols=vars(Site), rows=vars(variable), scales="free", space="free") +             
  scale_fill_manual("Taxon", values=c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', "black"),     
                    breaks = c("Araneae", "Pscoptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Acari", "Coleoptera", "Diptera"),  # Replace with your actual labels
                    labels = c("Araneae", "Psocoptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Acari", "Coleoptera", "Diptera"))+ # Replace with your desired legend labels
  scale_y_continuous(name="proportional zOTU representation", expand = c(0,0.1), breaks=seq(0,.60,.20))+
  #geom_text(aes(label = round(Arealog,1)),vjust=-.25, size=8) +
  #scale_x_discrete(breaks=rep$my_site, labels=rep$Arealog)+                 
  KipukaTheme +
  theme(axis.text.x = element_blank(), #element_text(angle=45, size=25, vjust=-.1, hjust=1),
        axis.text.y = element_text(size=25, vjust=0.5, hjust=0.5),
        axis.title.x=element_text(angle=0, size=30),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size=30), 
        legend.position = "top", 
        legend.box = "horizontal")

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
        strip.text.x = element_text(size=30),
        legend.position = "none"  )
                   
jpeg("Figures/Order_Representation.jpg", width=1500, height=2000)
plot_grid(a, b, nrow=2, rel_heights=c(1, .95))                   
dev.off()
