
#Set up my working environment
setwd("H:/My Drive/Kipuka/Data")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
library(data.table)
library(scales)
library(extrafont)

#Establish some color schemes up top to apply to all
#Colors are from color-blind friendly, rcartocolor "Safe" palette
SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainback"="#999933")
#Establish some themes up top to apply to all
KipukaTheme <- theme(axis.title=element_text(size=70), 
        axis.text = element_text(size=70, angle=45), 
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

richness <- read.csv("NatNonNat.csv")               
OTU <- richness[2:nrow(richness),29:814] 
#Make first column the row names
OTU["16",1]<-"ZotuID"
OTU["18",1]<-"3%OTU"
OTU <- OTU[-c(10),]
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)
#Filter out Araneae only!!!                
OTU<-OTU[OTU$Order=="Araneae",]               
OTU<-OTU[,-c(2:16)]
i <- c(2:ncol(OTU))                                  
OTU[i] <- apply(OTU[i], 2,            
                    function(x) as.numeric(as.character(x)))
sapply(OTU, class) 

#Now summarize native/non-native species per site using apply
#needs to be data.table for apply to work
OTU<-setDT(OTU)
otu<-OTU[, lapply(.SD, sum), by = INVNAT]
otu<- t(otu)
otu<-as.data.frame(otu)
names(otu) <- c("INV", "NAT")
otu<- otu[-1,]
i <- c(1:ncol(otu))                                  
otu[i] <- apply(otu[i], 2,            
                    function(x) as.numeric(as.character(x)))
otu$p_nat<-otu$NAT/(otu$INV+otu$NAT) #ranges from 42-100%
otu$p_non<-otu$INV/(otu$INV+otu$NAT) #ranges from 0-58%
otu$Site<-rownames(otu)
otu<-otu[,-c(1:2)]

#Convert wide to longform for plotting
otu<- melt(otu, id.vars=c("Site"))
otu <- otu[otu$Site != "1K08E" & otu$Site != "1K08C",]                      

#Join back other data needed to interpret sites
names(richness)<-richness[c(18),]
richness<-richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))                       

richness<-richness[1:11]
OTU<-merge(richness, otu, by.x="ID", by.y="Site")

#Order as I want panels to appear in plot
OTU$variable <- factor(OTU$variable, levels = rev(c("p_nat", "p_non")))                 
OTU$Area<-round(OTU$Area,10)    
plotA<-OTU                
plotA$Site<-gsub("Stainbeck","Stainback",as.character(plotA$Site)) 
plotA$Site <- factor(plotA$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava")))
                                           
####################################################################################################
#Also do in terms of species richness of nat/non-nat

richness <- read.csv("NatNonNat.csv")
                
OTU <- richness[2:nrow(richness),29:814] 
#Make first column the row names
OTU["16",1]<-"ZotuID"
OTU["18",1]<-"3%OTU"
OTU <- OTU[-c(10),]
rownames(OTU) <- OTU[,1]
OTU <- OTU[,-1]
OTU<- t(OTU)
OTU<-as.data.frame(OTU)
#Filter out Araneae only!!!                
OTU<-OTU[OTU$Order=="Araneae",]                  
OTU<-OTU[,-c(2:16)]
i <- c(2:ncol(OTU))                                  
OTU[i] <- apply(OTU[i], 2,            
                    function(x) as.numeric(as.character(x)))
sapply(OTU, class) 

#Now summarize native/non-native species per site using apply
#needs to be data.table for apply to work

for (I in 2:ncol(OTU)){
        for (J in 1:nrow(OTU)){
                if (OTU[J,I] > 0) {
                        OTU[J,I]<-1
                }   
                else if(OTU[J,I]==0){
                        OTU[J,I]<-0
               }
               else {
                       print("Something's wrong")
              }               
        }                
}        

OTU<-setDT(OTU)            
otu<-OTU[, lapply(.SD, sum), by = INVNAT]
otu<- t(otu)
otu<-as.data.frame(otu)
names(otu) <- c("INV", "NAT")
otu<- otu[-1,]
i <- c(1:ncol(otu))                                  
otu[i] <- apply(otu[i], 2,            
                    function(x) as.numeric(as.character(x)))
otu$p_nat<-otu$NAT/(otu$INV+otu$NAT)
otu$p_non<-otu$INV/(otu$INV+otu$NAT)
otu$Site<-rownames(otu)
otu<-otu[,-c(1:2)]

#Convert wide to longform for plotting
otu<- melt(otu, id.vars=c("Site"))
otu <- otu[otu$Site != "1K08E" & otu$Site != "1K08C",]                      

#Join back other data needed to interpret sites
names(richness)<-richness[c(18),]
richness<-richness[19:nrow(richness),]
richness$Area<-as.numeric(gsub(",","",as.character(richness$Area)))                       

richness<-richness[1:11]
OTU<-merge(richness, otu, by.x="ID", by.y="Site")

#Order as I want panels to appear in plot
OTU$Site <- factor(OTU$Site, levels = rev(c("Kona","Stainback",  "Center", "Edge", "Lava"))) 
OTU$variable <- factor(OTU$variable, levels = rev(c("p_nat", "p_non")))                                 
OTU$Area<-round(OTU$Area,10)                
plotB <-OTU[order(OTU$Site, OTU$Area),]
#Reindex data frame so that it plots this way
rownames(plotB) <- seq(1,nrow(plotB),1)                            

subsetB<-plotB[plotB$variable=="p_non",]
subsetB<-subsetB[subsetB$Site=="Center" | subsetB$Site== "Edge",] 
subsetB$Site<-gsub("Stainbeck","Stainback",as.character(subsetB$Site))               

###########################################################################################################
# TEST:  ANOVA to check whether proportion of non-native 3% OTUs ~ area type  ( lava, edge, center, Stainback, Kona)

# First, test assumptions:  Assumes normal distribution of data and equal variances between groups.               
# Check normality of residuals
test<-plotA[plotA$variable=="p_non",]                
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
        cat("Kruskal-Wallis test results", "\n")
        print(kruskal_result)
              pairwise_result <- pairwise.wilcox.test(test$value, test$Site, p.adj = "bonferroni")
              cat("Pairwise Wilcoxon test results")
              print(pairwise_result)
}  else { # Conduct the ANOVA                   
        anova_result <- aov(value ~ Site, data = test)
        print(paste0("ANOVA test results"))
        print(summary(anova_result))  
              # Post hoc Tukey's HSD test
              tukey_result <- TukeyHSD(anova_result)
              cat("Tukey's HSD test results")
              print(tukey_result)
}
                
###########################################################################################################
# TEST:  Linear regression for size vs. richness
      
for (i in 1:length(c("Center", "Edge"))){  
        type<-c("Center", "Edge")[i]
        print(type)
        test<-subsetB[subsetB$variable=="p_non" & subsetB$Site==type,]
        # First, test assumptions:  
        # Fit linear regression model
        linear_model <- lm(value ~ log10(Area), data = test)
        # Check assumptions
        # 1. Residuals vs Fitted Values Plot
        plot(residuals(linear_model) ~ fitted(linear_model), main = "Residuals vs Fitted", xlab = "Fitted Values",ylab = "Residuals")
        abline(h = 0, col = "red", lty = 2)
        # 2. Normal Q-Q Plot
        qqnorm(residuals(linear_model))
        qqline(residuals(linear_model), col = "red")
        # 3. Scale-Location (Spread-Location) Plot
        plot(sqrt(abs(residuals(linear_model))) ~ fitted(linear_model), main = "Scale-Location Plot", xlab = "Fitted Values", ylab = "sqrt(|Residuals|)")
        abline(h = 0, col = "red", lty = 2)
        # 4. Residuals vs Leverage Plot (Cook's distance)
        plot(hatvalues(linear_model), cooks.distance(linear_model), main = "Residuals vs Leverage", xlab = "Leverage", ylab = "Cook's distance")
        abline(h = 4/length(test$value), col = "red", lty = 2)
        # Print summary of the linear model
        print(paste0("linear regression for ", type))
        print(summary(linear_model))             
}
            
                
###########################################################################################################
# Now plot this
                
bp <- ggplot() + 
  geom_boxplot(data=plotA[plotA$variable=="p_non",],aes(x=Site, y=value, fill=Site), color="black", size=2)+
  scale_fill_manual(values=SiteColors) +
  labs(title="A. ", x="", y="Proportion of non-native species") +
  KipukaTheme +
  theme(strip.text = element_text(size = 40), 
        axis.text = element_text(angle=45, size=70, hjust=1, vjust=0), 
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
        legend.text=element_text(size=70), 
        legend.title = element_blank(),
       legend.position = "none", 
        plot.margin = margin(0.2,1,0,1.35, "cm"),
       legend.key.height = unit(3, 'cm'), 
        legend.key.width = unit(1, 'cm'))
                
B1<-ggplot() + 
  geom_point(data=subsetB,aes(x=as.numeric(Area), y=value, colour=Site), size=6, shape=15)+
  scale_colour_manual(values=SiteColors, limits=c("Center", "Edge")) +
  geom_smooth(method='lm', data=subsetB, aes(x=Area, y=value, colour=Site, fill=Site), size=1, alpha=0.20)+
  scale_fill_manual(values=SiteColors, limits=c("Center", "Edge")) +
  scale_y_continuous(name="Proportion of non-native OTUs")+              
  labs(title="B.", x="Kipuka area ("~m^2~")") +
  guides(colour="none", fill="none")+
  KipukaTheme +              
scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))  +                                    
  theme(strip.text = element_text(size = 40), 
        axis.text = element_text(angle=45, size=70, hjust=1, vjust=0), 
        axis.title.y = element_text(size=70), 
        panel.grid.major = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size=1),   
      panel.grid.minor = element_line(
        rgb(105, 105, 105, maxColorValue = 255),
        linetype = "dotted", 
        size = 0.5), 
       axis.title.x=element_text(size=70, margin=margin(-20,0,0,0)), 
        plot.title=element_text(size=70),  
        plot.margin = margin(0.2,1,0,1.35, "cm"))
                     
jpeg("../Figures/NatNonNat_AraneaeScatter_V2.jpg", width=2400, height=1200)
plot_grid(bp, B1, ncol=2, rel_widths=c(.7, 1))
dev.off()   
