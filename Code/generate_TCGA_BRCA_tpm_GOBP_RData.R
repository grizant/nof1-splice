## Store TCGA BRCA TPM as RData for convenience and computational efficiency
## subset to only genes in filtered GO-BP 2015
## AG Schissler
## Created 25 Jul 2017

################################
## 1. Create environment and read in counts
setwd("~/Dropbox/Splice-n-of-1-pathways")

load("Data/brca_iso_paired_tpm.RData")

################################
## 3. Load gene set definitions to subset iso_data (makes computation more manageable)
gobp <- read.delim2("~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt")

gobp_genes <- unique(gobp$symbol)
length(gobp_genes) ## 11530 genes

brca_iso_gobp_data <- brca_iso[brca_iso$geneSymbol %in% gobp_genes,]
nrow(brca_iso_gobp_data) ## 36710 isoforms
nrow(brca_iso)

## free some memory
rm(brca_iso)

################################
## 4. Load gene set definitions to subset iso_data (makes computation more manageable)

## save data frame as an R objects
save(brca_iso_gobp_data, file = "~/Dropbox/Splice-n-of-1-pathways/Data/brca_iso_gobp_data.RData")
