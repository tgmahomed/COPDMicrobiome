library("qiime2R")
library("phyloseq")
library("ggplot2")
library("vegan")
library("tidyverse")
library("ape")
library("RColorBrewer")
library("viridis")
library("ggsci")

S25_1=import_biom("F:/Virtual Machines/Shared/PhD and other work/NoHuman_Kraken/S25.biom", parseFunction = parse_taxonomy_greengenes)
S15_1=import_biom("F:/Virtual Machines/Shared/PhD and other work/NoHuman_Kraken/S15.biom", parseFunction = parse_taxonomy_greengenes)
S14_1=import_biom("F:/Virtual Machines/Shared/PhD and other work/NoHuman_Kraken/S14.biom", parseFunction = parse_taxonomy_greengenes)
S11_1=import_biom("F:/Virtual Machines/Shared/PhD and other work/NoHuman_Kraken/S11.biom", parseFunction = parse_taxonomy_greengenes)
S9_1=import_biom("F:/Virtual Machines/Shared/PhD and other work/NoHuman_Kraken/S9.biom", parseFunction = parse_taxonomy_greengenes)
S8_1=import_biom("F:/Virtual Machines/Shared/PhD and other work/NoHuman_Kraken/S8.biom", parseFunction = parse_taxonomy_greengenes)
sample_names(S25_1)="TG-M25"
sample_names(S15_1)="TG-M15"
sample_names(S14_1)="TG-M14"
sample_names(S11_1)="TG-M11"
sample_names(S9_1)="TG-M9"
sample_names(S8_1)="TG-M8"
mergeBiom=merge_phyloseq(S8_1, S9_1, S11_1, S14_1, S15_1, S25_1)

#Import Metadata file
Meta1=read.delim(file ="C:/Users/Tanweer/Documents/FilesForR/sample_metadata_Bac_vir_2.txt" , header = TRUE, sep = "\t", row.names = 1)
Meta2=sample_data(Meta1)

Virseq=merge_phyloseq(mergeBiom, Meta2)

#Plot Phyla
Virseq_phyla_glom=tax_glom(Virseq, "Phylum")
plot_bar(Virseq, fill = "Phylum")+ geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position = "stack") +  coord_flip() + theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle = 360), legend.title =element_text(colour = "black", size = 12, face = "bold"), legend.text =element_text(colour = "black", size = 12, face = "bold"), axis.title=element_text(colour = "black", size = 12, face = "bold"), legend.position = "bottom")
#Plot Family
plot_bar(Virseq, fill = "Family")+ geom_bar(aes(color=Family, fill=Family), stat = "identity", position = "stack") + theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle = 360), legend.title =element_text(colour = "black", size = 6, face = "bold"), legend.text =element_text(colour = "black", size = 8, face = "bold"), axis.title=element_text(colour = "black", size = 12, face = "bold"), legend.position = "bottom", legend.box = "horizontal") + scale_color_viridis(discrete = TRUE, option = "D") + scale_fill_viridis(discrete = TRUE)+ coord_flip() 

#write files
tax_table(Virseq)[is.na(tax_table(Virseq))] <- "Unclassified"
abudancedf_withtaxa_taxatable=data.frame(tax_table(Virseq))
write.csv(abudancedf_withtaxa_taxatable, file = "C:/Users/Tanweer/OneDrive/Documents/PhD/Virseq18Nov_taxatable.csv")

abudancedf_withtaxa_otutable=data.frame(otu_table(Virseq))
write.csv(abudancedf_withtaxa_otutable, file = "C:/Users/Tanweer/OneDrive/Documents/PhD/Virseq_18Nov_otutable.csv")

#realtive abudance
Virseq_RA = transform_sample_counts(Virseq, function(x) x / sum(x))
tax_table(Virseq_RA)[is.na(tax_table(Virseq_RA))] <- "Unclassified"
abudancedf_withtaxa_taxatable_1=data.frame(tax_table(Virseq_RA))
write.csv(abudancedf_withtaxa_taxatable_1, file = "C:/Users/Tanweer/OneDrive/Documents/PhD/VirseqRA18Nov_taxatable.csv")

abudancedf_withtaxa_otutable_1=data.frame(otu_table(Virseq_RA))
write.csv(abudancedf_withtaxa_otutable_1, file = "C:/Users/Tanweer/OneDrive/Documents/PhD/VirseqRA_18Nov_otutable.csv")


#Relative abundance_family
Virseq_fam_glom=tax_glom(Virseq, "Family")
Virseq_RA_1 = transform_sample_counts(Virseq_fam_glom, function(x) x / sum(x))
tax_table(Virseq_RA_1)[is.na(tax_table(Virseq_RA_1))] <- "Unclassified"
abudancedf_withtaxa_taxatable_2=data.frame(tax_table(Virseq_RA_1))
write.csv(abudancedf_withtaxa_taxatable_2, file = "C:/Users/Tanweer/OneDrive/Documents/PhD/VirseqRA18Nov_taxatable_Family.csv")

abudancedf_withtaxa_otutable_2=data.frame(otu_table(Virseq_RA_1))
write.csv(abudancedf_withtaxa_otutable_2, file = "C:/Users/Tanweer/OneDrive/Documents/PhD/VirseqRA_18Nov_otutable_Family.csv")

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  