library("qiime2R")
library("phyloseq")
library("ggplot2")
library("vegan")
library("tidyverse")
library("ape")
library("RColorBrewer")
library("viridis")

#Importing QIIME2 files
Physeq=qza_to_phyloseq(features ="C:/Users/Tanweer/Documents/FilesForR/table_deblur.qza", taxonomy ="C:/Users/Tanweer/Documents/FilesForR/taxonomy_1.qza")
sample_names(Physeq)=paste(sample_names(Physeq), "QIIME2", sep = "_")

OTU1=read.csv(file ="C:/Users/Tanweer/Documents/FilesForR/IS-Pro_OTUtable_6.csv" , header = TRUE, sep = ";", row.names =1)
OTU2=as.matrix(OTU1)
OTU3=otu_table(OTU2, taxa_are_rows = TRUE)
TAX1=read.csv(file = 'C:/Users/Tanweer/Documents/FilesForR/Tax_table_22June.csv', header = TRUE, sep = ";", row.names =1)
TAX2=tax_table(as.matrix(TAX1))

ISPhyseq=merge_phyloseq(OTU3, TAX2)

Merge_Physeq=merge_phyloseq(ISPhyseq, Physeq)


Meta1=read.delim(file ="C:/Users/Tanweer/Documents/FilesForR/Metadata_forComp_6.txt" , header = TRUE, sep = "\t", row.names = 1)
Meta2=sample_data(Meta1)

Merge_Physeq2=merge_phyloseq(Merge_Physeq, Meta2)
Merge_physeq3=subset_samples(Merge_Physeq2, sample_names(Merge_Physeq2)!="TG-M25_QIIME2")

Merge_Physeq4=subset_taxa(Merge_physeq3, Kingdom!="Archaea" & Family!="mitochondria" & Class!="Chloroplast")
tax_table(Merge_Physeq4)[is.na(tax_table(Merge_Physeq4))] <- "Unclassified"
Merge_Physeq_NA=subset_taxa(Merge_Physeq4, Order=="Unclassified")
Merge_Physeq_Method=merge_samples(Merge_Physeq_NA, "Method")
Merge_Physeq_NA_RA=transform_sample_counts(Merge_Physeq_NA, function(x) x/sum(x))

plot_bar(Merge_Physeq_Method)+ geom_bar(aes(fill=Phylum), stat = "identity", position = "stack") + theme(legend.position = "bottom", legend.text = element_text(size = 10))

#Merge_Physeq_NA_Phylum=subset_taxa(Merge_Physeq4, Phylum=="Unclassified")
Merge_Physeq_NA_Class=subset_taxa(Merge_Physeq4, Class=="Unclassified")
Merge_Physeq_NA_Order=subset_taxa(Merge_Physeq4, Order=="Unclassified")
Merge_Physeq_NA_Family=subset_taxa(Merge_Physeq4, Family=="Unclassified")
Merge_Physeq_NA_Genus=subset_taxa(Merge_Physeq4, Genus=="Unclassified")
Merge_Physeq_NA_Species=subset_taxa(Merge_Physeq4, Species=="Unclassified")

#Merge_Physeq_NA_Phylum_Method=merge_samples(Merge_Physeq_NA_Phylum, "Method")
Merge_Physeq_NA_Class_Method=merge_samples(Merge_Physeq_NA_Class, "Method")
Merge_Physeq_NA_Order_Method=merge_samples(Merge_Physeq_NA_Order, "Method")
Merge_Physeq_NA_Family_Method=merge_samples(Merge_Physeq_NA_Family, "Method")
Merge_Physeq_NA_Genus_Method=merge_samples(Merge_Physeq_NA_Genus, "Method")
Merge_Physeq_NA_Species_Method=merge_samples(Merge_Physeq_NA_Species, "Method")

Merge_Physeq_NA_Class_Method_RA=transform_sample_counts(Merge_Physeq_NA_Class_Method, function(x) x/sum(x))
plot_bar(Merge_Physeq_NA_Class_Method)+ geom_bar(aes(fill=Phylum), stat = "identity", position = "stack") + theme(legend.position = "bottom", legend.text = element_text(size = 10))
plot_bar(Merge_Physeq_NA_Class_Method_RA)+ geom_bar(aes(fill=Phylum), stat = "identity", position = "stack") + theme(legend.position = "bottom", legend.text = element_text(size = 10))

plot_bar(Merge_Physeq_NA_Class_Method)+ geom_bar(aes(fill=Class), stat = "identity", position = "stack") + theme(legend.position = "bottom", legend.text = element_text(size = 10))
plot_bar(Merge_Physeq_NA_Order_Method)+ geom_bar(aes(fill=Order), stat = "identity", position = "stack") + theme(legend.position = "bottom", legend.text = element_text(size = 10))
plot_bar(Merge_Physeq_NA_Family_Method)+ geom_bar(aes(fill=Family), stat = "identity", position = "stack") + theme(legend.position = "bottom", legend.text = element_text(size = 10))
plot_bar(Merge_Physeq_NA_Genus_Method)+ geom_bar(aes(fill=Genus), stat = "identity", position = "stack") + theme(legend.position = "bottom", legend.text = element_text(size = 10))
plot_bar(Merge_Physeq_NA_Species_Method)+ geom_bar(aes(fill=Species), stat = "identity", position = "stack") + theme(legend.position = "bottom", legend.text = element_text(size = 10))

#Changing aesthetics
p=plot_bar(Merge_Physeq_NA_Class_Method_RA)+ geom_bar(aes(fill=Phylum), stat = "identity", position = "stack") + theme(legend.position = "right", legend.text = element_text(size = 10), axis.text.x = element_text (size=10, face="bold", angle = 360, color="black"), axis.text.y =  element_text (size=10, face="bold", color="black"), axis.title = element_text(size = 12, face="bold"), strip.text.x = element_text(size = 14, color = "black", face = "bold"), legend.title = element_blank()) + scale_color_viridis(discrete = TRUE, option = "D") + scale_fill_viridis(discrete = TRUE)
p$data$Sample=as.character(p$data$Sample)
p$data$Sample=factor(p$data$Sample, level=c("Targeted metagenomics", "IS-Pro"))
print(p)

#
Merge_Physeq_NA_Class_Method_RA_IS_Pro=subset_samples(Merge_Physeq_NA_Class_Method_RA, method="IS-Pro")
