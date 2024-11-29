#Set up my working environment
setwd("H:/My Drive/Kipuka/Data")
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(tidyverse)
library(scales)
library(mgcv)
#install.packages("remotes")
#remotes::install_github("Jtrachsel/funfuns")
library(funfuns)
library(BiocManager)
library(vegan)
#BiocManager::install("phyloseq")
library(phyloseq)
library(tidyr)
library(dplyr)
library(ape)


richness <- read.csv("merged_by_site_2.csv")
phylo <- richness[c(8:10,12:17),32:ncol(richness)]
row.names(phylo)<-phylo$ZOTU
phylo<-phylo[-1]
phylo <- as.data.frame(t(phylo))
phylo <- phylo[, rev(seq_len(ncol(phylo)))]
# Reverse only the first two columns
phylo <- phylo[, c(2, 1, 3:ncol(phylo))]
phylo <- phylo %>% select(2, 1, everything()[-(1:2)])

# Create a hierarchical string to represent relationships
phylo$taxonomy <- apply(phylo[, -1], 1, paste, collapse = "/")
phylo$taxonomy <- as.factor(phylo$taxonomy)
# Extract unique taxonomies for tree construction
unique_taxonomies <- unique(phylo$taxonomy)
# Convert taxonomy strings into a tree using as.phylo
taxonomy_tree <- as.phylo(~taxonomy, data = phylo)
# Export tree to Newick format
write.tree(taxonomy_tree, file = "taxonomy_tree.newick")

# Optionally, visualize the tree
plot(taxonomy_tree)
