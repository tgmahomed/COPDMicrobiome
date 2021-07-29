library("qiime2R")
library("phyloseq")
library("ggplot2")
library("vegan")
library("tidyverse")
library("ape")
library("RColorBrewer")
library("viridis")


physeq=qza_to_phyloseq(features ="C:/Users/Tanweer/Documents/FilesForR/table_deblur.qza", tree ="C:/Users/Tanweer/Documents/FilesForR/unrooted-tree.qza", taxonomy ="C:/Users/Tanweer/Documents/FilesForR/taxonomy_1.qza", metadata ="C:/Users/Tanweer/Documents//FilesForR/sample-metadata_1.txt" )
physeq2=subset_taxa(physeq, Kingdom!="Archaea" & Family!="mitochondria" & Class!="Chloroplast")


physeq.distWUF=distance(physeq2, method="wunifrac")
physeq.distWUF.ord=ordinate(physeq2, method = "PCoA", distance = physeq.distWUF)
plot_ordination(physeq2, physeq.distWUF.ord, color="DiseaseState") + geom_point()+ stat_ellipse()  +   theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size = 10, face = "bold.italic"), axis.text.x = element_text (size=10, face="bold", angle = 360), axis.text.y =  element_text (size=10, face="bold.italic"), axis.title = element_text(size = 12, face="bold"), strip.text.x = element_text(size = 14, color = "black", face = "bold"))

plot_ordination(physeq2, physeq.distWUF.ord, color="DiseaseState") + geom_point()+ stat_ellipse(type='t',size =1)  +   theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size = 10, face = "bold.italic"), axis.text.x = element_text (size=10, face="bold", angle = 360), axis.text.y =  element_text (size=10, face="bold.italic"), axis.title = element_text(size = 12, face="bold"), strip.text.x = element_text(size = 14, color = "black", face = "bold")) 


