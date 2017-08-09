## Store TCGA BLCA TPM as RData for convenience and computational efficiency
## subset to only genes in filtered KEGG-BP 2015
## AG Schissler
## Created 25 Jul 2017

################################
## 1. Create environment and read in counts
setwd("~/Dropbox/Splice-n-of-1-pathways")

load("Data/blca_iso_paired_tpm.RData")

################################
## 3. Load gene set definitions to subset iso_data (makes computation more manageable)
kegg <- read.delim2("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")

kegg_genes <- unique(kegg$symbol)
length(kegg_genes) ## 5879 genes

blca_iso_kegg_data <- blca_iso[blca_iso$geneSymbol %in% kegg_genes,]
nrow(blca_iso_kegg_data)
kegg_measured <- length(unique(blca_iso_kegg_data$geneSymbol))
nrow(blca_iso)

## free some memory
rm(blca_iso)

## how many genes isoforms where filtered?
iso_range <- c(2,30)

iso_list <- split(blca_iso_kegg_data$isoform_id, factor(blca_iso_kegg_data$geneSymbol))
str(tmp_gene <- iso_list[[1]])
to_keep <- unlist(lapply(iso_list, function(tmp_gene) {
    tmp_num <- length(tmp_gene)
    tmp_num >= iso_range[1] & tmp_num <= iso_range[2]
}))
kegg_measured - sum(!to_keep) ## filter 1624 genes
kegg_genes <- sum(!to_keep)
nrow(blca_iso_kegg_data[blca_iso_kegg_data$geneSymbol %in% names(which(to_keep)),]) ## 17088

################################
## 4. Load gene set definitions to subset iso_data (makes computation more manageable)

## save data frame as an R objects
save(blca_iso_kegg_data, file = "~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_kegg_data.RData")
