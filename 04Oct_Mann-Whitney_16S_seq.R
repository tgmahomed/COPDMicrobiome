library("qiime2R")
library("phyloseq")
library("ggplot2")
library("vegan")
library("tidyverse")
library("ape")
library("RColorBrewer")
library("viridis")

#Importing QIIME2 files
physeq=qza_to_phyloseq(features ="C:/Users/Tanweer/Documents/FilesForR/table_deblur.qza", tree ="C:/Users/Tanweer/Documents/FilesForR/unrooted-tree.qza", taxonomy ="C:/Users/Tanweer/Documents/FilesForR/taxonomy_1.qza", metadata ="C:/Users/Tanweer/Documents/FilesForR/sample-metadata_1.txt" )
physeq2=subset_taxa(physeq, Kingdom!="Archaea" & Family!="mitochondria" & Class!="Chloroplast")


#Calculate Richness
richness=estimate_richness(physeq2)
write.table(richness, "C:/Users/Tanweer/Documents/FilesForR/AlphaDiveristy_physeq2.txt", sep = "\t")
#Merge phyloseq for method and calculate richness
#Calcualte Mann-Whitney
pairwise.wilcox.test(richness$Chao1, sample_data(physeq2)$DiseaseState)
pairwise.wilcox.test(richness$Simpson, sample_data(physeq2)$DiseaseState)
