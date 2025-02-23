# Set up my working environment
setwd("H:/My Drive/Kipuka/Data")
library(reshape2)
library(tidyverse)
library(vegan)
library(dplyr)

otu <- read.csv("OTUs.csv")
richness <- read.csv("merged_by_site_2.csv")

# Grab site name, site type, and coordinate data for later use
site_data<-as.data.frame(richness[19:nrow(richness), c(9,27,28,29)])
names(site_data)<-c("Site", "my_ID", "Latitude", "Longitude")
site_data$Site<-gsub("Stainbeck","Stainback",as.character(site_data$Site))

# Clean up OTU table for use in distance matrix calculation
OTU <- otu[14:nrow(otu), 29:ncol(otu)] 
rownames(OTU) <- OTU[,1]
OTU <- as.data.frame(t(OTU[,-1]))
names(OTU)[1:2] <- c("OTU", "zOTU")
# Exclude blank (NA or empty string) values 
OTUtoKeep <- OTU[grepl("OTU", OTU[[1]]), ]
OTUtoKeep[3:ncol(OTUtoKeep)] <- lapply(OTUtoKeep[3:ncol(OTUtoKeep)], as.numeric)

# OTU
OTU3 <- OTUtoKeep %>%
  group_by(OTU) %>%
  summarise(across(2:(ncol(OTUtoKeep)-1), sum, na.rm = TRUE)) # Sum columns 3 to last
OTU3 <- OTU3[,-1]
OTU3 <- as.data.frame(t(OTU3))
OTU3[] <- lapply(OTU3, as.numeric)     

# Create B-C distance matrix
OTUbeta <- vegdist(OTU3, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
# Now, convert to long-form, and export for use in questions of beta diversity. 
OTUbeta<-as.matrix(OTUbeta)
OTUbeta<-as.data.frame(OTUbeta)
dist_long_3OTU <- melt(as.matrix(OTUbeta))
write.csv(dist_long_3OTU, "OTU3_Bray.csv", quote=F, row.names=F)

# Create Jaccard distance matrix
OTUbeta <- vegdist(OTU3, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
# Now, convert to long-form, and export for use in questions of beta diversity. 
OTUbeta<-as.matrix(OTUbeta)
OTUbeta<-as.data.frame(OTUbeta)
dist_long_3OTU <- melt(as.matrix(OTUbeta))
write.csv(dist_long_3OTU, "OTU3_jaccard.csv", quote=F, row.names=F)

# Mantel test for spatial autocorrelation based on 3% radius OTU matrix 
# FIrst, join the OTU3 data with the site data... 
mantel<-OTU3
mantel$my_ID<-rownames(OTU3)
mantel<-merge(mantel, site_data, by="my_ID")
mantel<-mantel[,-1]
# Now, iterate through each site type
for (X in 1:length(unique(mantel$Site))){
  SITE <- unique(mantel$Site)[X]
  # Subset OTU matrix and coordinates for that site type
  INDICES <- which(mantel$Site==SITE)
  SUB_MAT <- mantel[INDICES,]
  SUB_MAT[1:ncol(SUB_MAT)] <- lapply(SUB_MAT[1:ncol(SUB_MAT)], as.numeric)
  SUB_COORDS <- mantel[INDICES, c("Latitude", "Longitude")]
  coords <- cbind(as.numeric(SUB_COORDS$Longitude), as.numeric(SUB_COORDS$Latitude))
  # Create a distance matrix from subset coordinates
  dist_matrix <- dist(coords)  
  # Create a B-C distance matrix from OTUs
  OTUbeta <- vegdist(SUB_MAT, method="bray", binary=FALSE, na.rm=T)
  print(paste0("3% OTU and Bray-Curtis at site type ", SITE))
  print(mantel(OTUbeta, dist_matrix))
}




# zOTU
# Identify the values in column 1 that occur more than once
#values_to_keep <- names(which(table(OTUtoKeep[[1]]) > 1))

# Subset the dataframe to keep only rows where column 1 matches those values
OTUtoKeep_filtered <- OTUtoKeep#[OTUtoKeep[[1]] %in% values_to_keep, ]
OTUtoKeep_filtered <- OTUtoKeep_filtered[,-c(1,2)]
OTUtoKeep_filtered <- as.data.frame(t(OTUtoKeep_filtered))
OTUtoKeep_filtered[] <- lapply(OTUtoKeep_filtered, as.numeric)

zOTUbeta <- vegdist(OTUtoKeep_filtered, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-as.data.frame(zOTUbeta)
dist_long <- melt(as.matrix(zOTUbeta))
write.csv(dist_long, "zOTU_Bray.csv", quote=F, row.names=F)

zOTUbeta <- vegdist(OTUtoKeep_filtered, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-as.data.frame(zOTUbeta)
dist_long_zOTU <- melt(as.matrix(zOTUbeta))
write.csv(dist_long, "zOTU_jaccard.csv", quote=F, row.names=F)

# Mantel test for spatial autocorrelation based on zOTU matrix 
# FIrst, join the OTU3 data with the site data... 
zantel<-OTUtoKeep_filtered
zantel$my_ID<-rownames(OTUtoKeep_filtered)
zantel<-merge(zantel, site_data, by="my_ID")
zantel<-zantel[,-1]
# Now, iterate through each site type
for (X in 1:length(unique(zantel$Site))){
  SITE <- unique(zantel$Site)[X]
  # Subset OTU matrix and coordinates for that site type
  INDICES <- which(zantel$Site==SITE)
  SUB_MAT <- zantel[INDICES,]
  SUB_MAT[1:ncol(SUB_MAT)] <- lapply(SUB_MAT[1:ncol(SUB_MAT)], as.numeric)
  SUB_COORDS <- zantel[INDICES, c("Latitude", "Longitude")]
  coords <- cbind(as.numeric(SUB_COORDS$Longitude), as.numeric(SUB_COORDS$Latitude))
  # Create a distance matrix from subset coordinates
  dist_matrix <- dist(coords)  
  # Subset 3OTU matrix also for that site type
  INDICES <- which(mantel$Site==SITE)
  MAT3 <- mantel[INDICES,]
  MAT3[1:ncol(MAT3)] <- lapply(MAT3[1:ncol(MAT3)], as.numeric)
  # Create a Jaccard distance matrix from OTUs
  zOTUdist <- vegdist(SUB_MAT, method="jaccard", binary=FALSE, na.rm=T)
  OTUdist <- vegdist(MAT3, method="jaccard", binary=FALSE, na.rm=T)
  print(paste0("zOTU and Jaccard at site type ", SITE))
  print(mantel(dist_matrix, zOTUdist))
  partial_mantel_result <- mantel.partial(zOTUdist, dist_matrix, OTUdist, method = "pearson", permutations = 999)
  print(partial_mantel_result)
}



# Now, we need a 3% OTU Bray-Curtis per-taxon... 
richness <- read.csv("merged_by_site_2.csv")
OTUtoKeep<-as.data.frame(t(richness[c(12,17:nrow(richness)), 33:ncol(richness)]))
names(OTUtoKeep)[1:3]<- c("Order","OTU", "zOTU")
names(OTUtoKeep)[4:ncol(OTUtoKeep)] <- richness[19:nrow(richness),1]
OTUtoKeep <- OTUtoKeep[OTUtoKeep$Order %in% c("Araneae", "Lepidoptera", "Coleoptera", "Diptera", "Psocoptera", "Hemiptera"),]
# Exclude blank (NA or empty string) values 
OTUtoKeep <- OTUtoKeep[grepl("OTU", OTUtoKeep[[2]]), ]
OTUtoKeep[4:ncol(OTUtoKeep)] <- lapply(OTUtoKeep[4:ncol(OTUtoKeep)], as.numeric)

all_dist_long <- list()

for (i in 1:length(unique(OTUtoKeep$Order))){
ORDER <- unique(OTUtoKeep$Order)[i]
OTU3 <- OTUtoKeep[OTUtoKeep$Order==ORDER,] %>%
  group_by(OTU) %>%
  summarise(across(3:(ncol(OTUtoKeep)-1), sum, na.rm = TRUE)) # Sum columns 3 to last
OTU3 <- OTU3[,-1]
OTU3 <- as.data.frame(t(OTU3))
OTU3[] <- lapply(OTU3, as.numeric)
zOTUbeta <- vegdist(OTU3, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-as.data.frame(zOTUbeta)
# Melt the matrix and rename the value column
  # Melt the distance matrix and rename the value column
  dist_long <- melt(as.matrix(zOTUbeta))
  value_column_name <- paste0(ORDER, "_value")
  colnames(dist_long)[3] <- value_column_name
  
  # Store the data frame in a list
  all_dist_long[[i]] <- dist_long
}

# Merge all data frames by Var1 and Var2
merged_dist_long <- Reduce(function(x, y) full_join(x, y, by = c("Var1", "Var2")), all_dist_long)
write.csv(merged_dist_long, "OTU3_Bray_Order.csv", quote=F, row.names=F)



                           

# Now, we need a zOTU Jaccard per-taxon... 
                           
richness <- read.csv("merged_by_site_2.csv")
OTUtoKeep<-as.data.frame(t(richness[c(12,17:nrow(richness)), 33:ncol(richness)]))
names(OTUtoKeep)[1:3]<- c("Order","OTU", "zOTU")
names(OTUtoKeep)[4:ncol(OTUtoKeep)] <- richness[19:nrow(richness),1]
OTUtoKeep <- OTUtoKeep[OTUtoKeep$Order %in% c("Araneae", "Lepidoptera", "Coleoptera", "Diptera", "Psocoptera", "Hemiptera"),]
# Identify the values in OTU column that occur more than once
values_to_keep <- names(which(table(OTUtoKeep[[2]]) > 1))                                  
# Exclude blank (NA or empty string) values 
OTUtoKeep <- OTUtoKeep[grepl("OTU", OTUtoKeep[[2]]), ]
OTUtoKeep[4:ncol(OTUtoKeep)] <- lapply(OTUtoKeep[4:ncol(OTUtoKeep)], as.numeric)
OTUtoKeep_filtered <- OTUtoKeep[OTUtoKeep[[2]] %in% values_to_keep, ]   
OTUtoKeep_filtered <- OTUtoKeep_filtered[,-c(2,3)]                                             

all_dist_long <- list()

for (i in 1:length(unique(OTUtoKeep_filtered$Order))){
ORDER <- unique(OTUtoKeep_filtered$Order)[i]
zOTU <- as.data.frame(OTUtoKeep_filtered[OTUtoKeep_filtered$Order==ORDER,])
zOTU <- zOTU[,-1]
zOTU <- as.data.frame(t(zOTU))
zOTU[] <- lapply(zOTU, as.numeric)
zOTUbeta <- vegdist(zOTU, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
zOTUbeta<-as.matrix(zOTUbeta)
zOTUbeta<-as.data.frame(zOTUbeta)
# Melt the matrix and rename the value column
  # Melt the distance matrix and rename the value column
  dist_long <- melt(as.matrix(zOTUbeta))
  value_column_name <- paste0(ORDER, "_value")
  colnames(dist_long)[3] <- value_column_name  
  # Store the data frame in a list
  all_dist_long[[i]] <- dist_long
}

# Merge all data frames by Var1 and Var2
merged_dist_long <- Reduce(function(x, y) full_join(x, y, by = c("Var1", "Var2")), all_dist_long)
write.csv(merged_dist_long, "zOTU_Jaccard_Order.csv", quote=F, row.names=F)
