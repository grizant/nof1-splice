## Store TCGA COAD TPM as RData for convenience and computational efficiency
## subset to only genes in filtered KEGG-BP 2015
## AG Schissler
## Created 25 Jul 2017

################################
## 1. Create environment and read in counts
setwd("~/Dropbox/Splice-n-of-1-pathways")

load("Data/coad_iso_paired_tpm.RData")

################################
## 3. Load gene set definitions to subset iso_data (makes computation more manageable)
kegg <- read.delim2("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")

kegg_genes <- unique(kegg$symbol)
length(kegg_genes) ## 5879 genes

coad_iso_kegg_data <- coad_iso[coad_iso$geneSymbol %in% kegg_genes,]
nrow(coad_iso_kegg_data)
nrow(coad_iso)

## free some memory
rm(coad_iso)

################################
## 4. Load gene set definitions to subset iso_data (makes computation more manageable)

## save data frame as an R objects
save(coad_iso_kegg_data, file = "~/Dropbox/Splice-n-of-1-pathways/Data/coad_iso_kegg_data.RData")
