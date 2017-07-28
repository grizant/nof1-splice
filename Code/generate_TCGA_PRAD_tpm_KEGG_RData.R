## Store TCGA PRAD TPM as RData for convenience and computational efficiency
## subset to only genes in filtered KEGG-BP 2015
## AG Schissler
## Created 25 Jul 2017

################################
## 1. Create environment and read in counts
setwd("~/Dropbox/Splice-n-of-1-pathways")

load("Data/prad_iso_paired_tpm.RData")

length(unique(prad_iso$geneSymbol))

################################
## 3. Load gene set definitions to subset iso_data (makes computation more manageable)
kegg <- read.delim2("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")

kegg_genes <- unique(kegg$symbol)
length(kegg_genes) ## 5879 genes

prad_iso_kegg_data <- prad_iso[prad_iso$geneSymbol %in% kegg_genes,]
nrow(prad_iso_kegg_data)
nrow(prad_iso)

## free some memory
rm(prad_iso)

################################
## 4. Load gene set definitions to subset iso_data (makes computation more manageable)

## save data frame as an R objects
save(prad_iso_kegg_data, file = "~/Dropbox/Splice-n-of-1-pathways/Data/prad_iso_kegg_data.RData")
