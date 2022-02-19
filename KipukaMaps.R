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

#First, rainfall showing inset site polygons

rain<-raster("AET_mm_ann_hr_raster")
#rain <- readOGR("Mapping/Annual_Rainfall_(mm)/Annual_Rainfall_(mm).shp") 

ggplot() + 
  geom_polygon(data=rain, aes(long, lat))  #, group = group, fill = hole 



coords <- read.csv("merged_by_site_2.csv")
ggplot() + 
  #geom_polygon(data=rain, aes(long, lat)) + #, group = group, fill = hole 
  geom_point(data=coords, aes(lon, lat))
