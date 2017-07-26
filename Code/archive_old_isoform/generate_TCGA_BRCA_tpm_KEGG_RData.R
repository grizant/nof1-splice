## Store TCGA BRCA TPM as RData for convenience and computational efficiency
## subset to only genes in filtered KEGG-BP 2015
## AG Schissler
## Created 14 Apr 2017
## Last modified 14 Apr 2017

################################
## 1. Create environment and read in counts
setwd("~/Dropbox/Splice-n-of-1-pathways/UofUtah")

## read in FPKM Isoform data
library(data.table)
norm_data <- fread(file = "Total_TPM_normal.txt", sep = "\t", header = T, stringsAsFactors = F)
tumor_data <- fread(file = "Total_TPM_tumor.txt", sep = "\t", header = T, stringsAsFactors = F)

## check normalization
## colSums(norm_data[,-(1:4)])
## colSums(tumor_data[,-(1:4)])

iso_data <- merge(norm_data, tumor_data)

################################
## 2. Rename using T and N
names(iso_data) <- gsub("01[A-Z]","T", names(iso_data))
names(iso_data) <- gsub("11[A-Z]","N", names(iso_data))

################################
## 3. Load gene set definitions to subset iso_data (makes computation more manageable)
kegg <- read.delim2("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")

kegg_genes <- unique(kegg$symbol)
length(kegg_genes) ## 11530 genes

iso_kegg_data <- iso_data[iso_data$Gene_symbol %in% kegg_genes,]
nrow(iso_kegg_data)
nrow(iso_data)

## free some memory
rm(iso_data)

## remove the data.table class for convenience
class(iso_kegg_data) <- "data.frame"

################################
## 4. Load gene set definitions to subset iso_data (makes computation more manageable)

## save data frame as an R objects
save(iso_kegg_data, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_tpm_KEGG.RData")
