library("qiime2R")
library("phyloseq")
library("ggplot2")
library("vegan")
library("tidyverse")
library("ape")
library("RColorBrewer")
library("viridis")
library("DESeq2")

Physeq=qza_to_phyloseq(features="C:/Users/tgmah/OneDrive/Documents/PhD/R Files/New folder/table_deblur.qza", taxonomy ="C:/Users/tgmah/OneDrive/Documents/PhD/R Files/New folder/taxonomy_1.qza")
Meta1=read.delim(file ="C:/Users/tgmah/OneDrive/Documents/PhD/R Files/New folder/sample_metadata_Bac_vir_2.txt" , header = TRUE, sep = "\t", row.names = 1)
Meta2=sample_data(Meta1)
Bacseq=merge_phyloseq(Physeq, Meta2)

Di_des=phyloseq_to_deseq2(Bacseq, ~DiseaseState)
Di_des_1=DESeq(Di_des)


resultsNames(Di_des_1)
resdf=as.data.frame(DESeq2::results(Di_des_1, format = "DataFrame", name="DiseaseState_Stable_vs_Exacerbation"))
resdf_2=results(Di_des_1, contrast = c("DiseaseState", "Exacerbation", "Stable"))

res=results(Di_des_1, cooksCutoff = FALSE)
sigtab=res[which(res$padj <0.2), ]
sigtab_2=cbind(as(sigtab, "data.frame"), as(tax_table(Bacseq)[rownames(sigtab), ], "matrix"))
head(sigtab_2)
dim(sigtab_2)

head(sigtab [order(sigtab$log2FoldChange ), ] )


library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname= "Set1", ...) {scale_fill_brewer(palette = palname, ...)}

x=tapply(sigtab_2$log2FoldChange, sigtab_2$Phylum, function(x) max(x))
x=sort(x, TRUE)

x=tapply(sigtab_2$log2FoldChange, sigtab_2$Genus, function(x) max(x))
x=sort(x, TRUE)
sigtab_2$Genus=factor(as.character(sigtab_2$Genus), levels=names(x))
ggplot(sigtab_2, aes(x=Genus, y=log2FoldChange)) + geom_bar(stat = "identity", position=position_dodge(), aes(fill=Phylum)) + geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2, position = position_dodge(.9))+ theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))  + coord_flip() + geom_hline(yintercept = 0.0, color = "black", size = 1) +   theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size = 18, face = "bold.italic"), axis.text.x = element_text (size=18, face="bold", angle = 360), axis.text.y =  element_text (size=18, face="bold.italic"), axis.title = element_text(size = 22, face="bold"), strip.text.x = element_text(size = 24, color = "black", face = "bold"), plot.title = element_text(size = 28, color = "black", face = "bold", hjust = 0.5) 
) + scale_color_viridis(discrete = TRUE, option = "D") + scale_fill_viridis(discrete = TRUE)   + ggtitle("Stable vs Exacerbated COPD patients") 

