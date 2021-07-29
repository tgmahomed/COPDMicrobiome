#Install packages
.cran_packages <- c("tidyverse", "cowplot", "picante", "vegan", "HMP", "dendextend", "rms", "devtools")
.bioc_packages <- c("phyloseq", "DESeq2", "microbiome", "metagenomeSeq", "ALDEx2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(.bioc_packages, version = "3.9")
devtools::install_github("adw96/breakaway")
devtools::install_github(repo = "UVic-omics/selbal")


#Load packages
library("qiime2R")
library("phyloseq")
library("ggplot2")
library("vegan")
library("tidyverse")
library("ape")
library("RColorBrewer")
library("viridis")
library("microbiome")
library("picante")
library("ALDEx2")
library("HMP")
library("dendextend")
library("selbal")
library("rms")
library("breakaway")
library("dplyr")
library("tibble")

#Load phyloseq object
Physeq=qza_to_phyloseq(features="C:/Users/tgmah/OneDrive/Documents/PhD/R Files/New folder/table_deblur.qza", taxonomy ="C:/Users/tgmah/OneDrive/Documents/PhD/R Files/New folder/taxonomy_1.qza")
Meta1=read.delim(file ="C:/Users/tgmah/OneDrive/Documents/PhD/R Files/New folder/sample_metadata_Bac_vir_2.txt" , header = TRUE, sep = "\t", row.names = 1)
Meta2=sample_data(Meta1)
Bacseq=merge_phyloseq(Physeq, Meta2)


#Attempt1
aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(Bacseq)), phyloseq::sample_data(Bacseq)$DiseaseState, test="t", effect = TRUE, denom="iqlr")


ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)


#
taxa_info <- data.frame(tax_table(Bacseq))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
