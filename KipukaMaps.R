#Project: Kipukas               #
#Script: KipukaMaps.R     #
#Author: Emma Steigerwald       #
#Date:18 Feb 2022               #
#################################

#Set up my working environment
setwd("G:/My Drive/Kipuka")
library(ggplot2)
library(extrafont)
library(sp)
library(rgdal)
library(raster)
font_import()
library(viridis)

SiteColors <- c("Center" = "#332288", "Edge" = "#6699CC", "Lava"="#888888", "Kona"="#117733", "Stainbeck"="#999933")
#Establish some themes up top to apply to all
KipukaTheme <- theme(axis.title=element_text(size=30), 
        axis.text = element_text(size=25), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        legend.text=element_text(size=25), 
        legend.key.height = unit(1, "cm"), 
        legend.key.width = unit(1.5,"cm"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(),
        legend.title = element_text(size=25), 
        text = element_text(family = "serif")) 



#First, rainfall showing inset site polygons

rain<-raster("Mapping/HawaiiRFGrids_mm/rf_mm_bi_ann/w001001.adf")
rain_spdf <- as(rain, "SpatialPixelsDataFrame")
rain_df <- as.data.frame(rain_spdf)
colnames(rain_df) <- c("value", "x", "y")

kona<-readOGR("Mapping/extent/Kona_extent.shp")
stainbeck<-readOGR("Mapping/extent/Stainbeck_extent.shp")
kipuka<-readOGR("Mapping/extent/Kipuka_extent.shp")

jpeg("Figures/BigIsland.jpg", width=1000, height=1000)
ggplot() +  
  geom_tile(data=rain_df, aes(x=x, y=y, fill=value^.5)) +  
  scale_fill_viridis_c(direction=-1, breaks=c(2000^.5, 4000^.5, 6000^.5), labels=c(2000, 4000, 6000)) +
  geom_polygon(data=kona, color="black", lwd=2, color="black", fill=NA, mapping=aes(long, lat))+
  geom_polygon(data=kipuka, color="black", lwd=2, color="black", fill=NA,  mapping=aes(long, lat))+
  geom_polygon(data=stainbeck, color="black", lwd=2, color="black", fill=NA,  mapping=aes(long, lat))+  
  KipukaTheme
dev.off()


coords <- read.csv("merged_by_site_2.csv")
ggplot() + 
  #geom_polygon(data=rain, aes(long, lat)) + #, group = group, fill = hole 
  geom_point(data=coords, aes(lon, lat))
