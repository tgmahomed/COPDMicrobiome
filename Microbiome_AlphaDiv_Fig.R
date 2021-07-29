#Load packages
library("qiime2R")
library("phyloseq")
library("ggplot2")
library("vegan")
library("tidyverse")
library("ape")
library("RColorBrewer")
library("viridis")
library("ggsci")


#import bacteriome data
Physeq=qza_to_phyloseq(features="C:/Users/Tanweer/Documents/FilesForR/table_deblur.qza", taxonomy ="C:/Users/Tanweer/Documents/FilesForR/taxonomy_1.qza")
Meta1=read.delim(file ="C:/Users/Tanweer/Documents/FilesForR/sample_metadata_Bac_vir_3.txt" , header = TRUE, sep = "\t", row.names = 1)
Meta2=sample_data(Meta1)
Bacseq=merge_phyloseq(Physeq, Meta2)


plot_richness(Bacseq, x= "DiseaseState", color = "DiseaseState", measures = c("Chao1", "Simpson")) + geom_boxplot() + scale_color_brewer(palette = "Paired") +   theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size = 10, face = "bold.italic"), axis.text.x = element_text (size=10, face="bold", angle = 360), axis.text.y =  element_text (size=10, face="bold.italic"), axis.title = element_text(size = 12, face="bold"), strip.text.x = element_text(size = 14, color = "black", face = "bold"))

p <- plot_richness(Bacseq, "DiseaseState", "Participant", measures = c("Chao1", "Simpson"))
(p <- p + geom_boxplot(data = p$data, aes(x = DiseaseState, y = value, color = NULL),  alpha = 0.1)) +   theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size = 10, face = "bold.italic"), axis.text.x = element_text (size=10, face="bold", angle = 360), axis.text.y =  element_text (size=10, face="bold.italic"), axis.title = element_text(size = 12, face="bold"), strip.text.x = element_text(size = 14, color = "black", face = "bold"))  

richness=estimate_richness(Bacseq)
Bacseq_Stable=subset_samples(Bacseq, DiseaseState=="Stable")
Bacseq_Ex=subset_samples(Bacseq, DiseaseState=="Exacerbation")

mean(estimate_richness(Bacseq_Stable, measures = c("Chao1"))[,1])
mean(estimate_richness(Bacseq_Ex, measures = c("Chao1"))[,1])
mean(estimate_richness(Bacseq_Stable, measures = c("Simpson"))[,1])
mean(estimate_richness(Bacseq_Ex, measures = c("Simpson"))[,1])

median(estimate_richness(Bacseq_Stable, measures = c("Chao1"))[,1])
median(estimate_richness(Bacseq_Ex, measures = c("Chao1"))[,1])
median(estimate_richness(Bacseq_Stable, measures = c("Simpson"))[,1])
median(estimate_richness(Bacseq_Ex, measures = c("Simpson"))[,1])

IQR(estimate_richness(Bacseq_Stable, measures = c("Chao1"))[,1])
IQR(estimate_richness(Bacseq_Ex, measures = c("Chao1"))[,1])
IQR(estimate_richness(Bacseq_Stable, measures = c("Simpson"))[,1])
IQR(estimate_richness(Bacseq_Ex, measures = c("Simpson"))[,1])
