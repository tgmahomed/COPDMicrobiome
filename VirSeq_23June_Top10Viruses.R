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
Virseq_RA = transform_sample_counts(Virseq, function(x) x / sum(x))
Virseq_Family_glom=tax_glom(Virseq, "Family")
Virseq_Family_RA=transform_sample_counts(Virseq_Family_glom, function(x) x / sum(x))


df <- psmelt(Virseq_Family_RA)

top_family <- df %>% group_by(Sample, Family) %>% summarize(Mean = mean(Abundance)) %>%  arrange(-Mean)


top10 <- top_family$Family[1:10]
df0 <- df %>%
  mutate(Family= fct_other(Family, top10))
ggplot(df0, aes(Sample, Abundance, fill = Family)) +
  geom_col(aes(color=Family, fill=Family), stat = "identity", position = "stack") + theme(axis.text.y = element_text(colour = "black", size = 20, face = "bold"), axis.text.x = element_text(colour = "black", size = 20, face = "bold", angle = 360), legend.title =element_text(colour = "black", size = 20, face = "bold"), legend.text =element_text(colour = "black", size = 20, face = "bold" ), axis.title=element_text(colour = "black", size = 20, face = "bold"), legend.position = "bottom", legend.box = "horizontal") + coord_flip() 

