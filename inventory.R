# Summarizing data by counts of zOTUs and OTUs per taxonomic group considered #
# Author: Emma Steigerwald                                                    #
###############################################################################

setwd("G:/My Drive/Kipuka/Data")
library(stringr)
library(dplyr)

inv<-read.csv("inventory.csv", header=F)
inv<-t(inv)
colnames(inv) <- as.character(inv[1, ])  # Extract the values from the first row
inv <- inv[-1, , drop = FALSE]
inv<-as.data.frame(inv)
names(inv)[1:2]<-c("Class1", "Class2")
inv <- inv[!is.na(inv[, 1]), , drop = FALSE]
names(inv)[8]<-"threeOTU"
inv$fullName<-as.character(paste0(inv$Class1, "_", inv$Class2, "_", inv$Order, "_", inv$Family, "_", inv$Genus, "_", inv$Species))
inv_split <- data.frame(do.call(rbind, str_split(unique(inv$fullName), "_"))) # Split the fullName column into separate columns using the character '_'
# Combine the original data frame and the new columns
unique <- cbind(unique(inv$fullName), inv_split)
colnames(unique) <- c("fullName", "Class1", "Class2", "Order", "Family", "Genus", "Species") # Rename the new columns

unique <- unique %>%
  mutate(zOTUcounts = NA, threeCounts = NA)

for (i in 1:nrow(unique)) { # Iterate over each row in unique
  OR <- unique$Order[i]
  # Sum how many unique values of inv$ID occur for the specific fullName
  total_ID <- length(unique(inv[inv$Order == OR, "ID"]))
  # Sum how many unique values of inv$threeOTU occur for the specific fullName
  total_threeOTU <- length(unique(inv[inv$Order == OR, "threeOTU"]))
  # Assign the calculated values to the corresponding rows in unique
  unique[i, "zOTUcounts"] <- total_ID
  unique[i, "threeCounts"] <- total_threeOTU
}

# Check that nuber of zOTUs is always greater than or equal to number 3% radius OTUs
all(unique$zOTUcounts >= unique$threeCounts)

write.csv(unique, "inventory_produced.csv", quote=F)
