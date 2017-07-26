## Explore KEGG for validation
## AG Schissler
## Created 18 Jul 2017

############################################################
## 1. read in gene sets

library(nof1)

annot_file <- "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt"
desc_file <- "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt"

annot_data <- nof1::read_gene_set(annot_file)
desc_data <- nof1::read_gene_set(desc_file, quote = "")

desc_data[grep("cancer", desc_data$description),]
