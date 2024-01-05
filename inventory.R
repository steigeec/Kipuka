# Summarizing data by counts of zOTUs and OTUs per taxonomic group considered #
# Author: Emma Steigerwald                                                    #
###############################################################################

inv<-read.csv("inventory.csv")
inv<-t(inv)
rownames(inv)<-inv[[1]]
inv<-inv[-1]

# Create a blank dataframe for summary info
bio <- data.frame(
  Class1 = factor(rep(NA, length(unique_species))),
  Class2 = factor(rep(NA, length(unique_species))),
  Order = factor(rep(NA, length(unique_species))),
  Family = factor(rep(NA, length(unique_species))),
  Genus = factor(rep(NA, length(unique_species))),
  Species = factor(unique(inv$species)),
  zOTUcounts = factor(rep(NA, length(unique_species))),
  threeCounts = factor(rep(NA, length(unique_species))))

