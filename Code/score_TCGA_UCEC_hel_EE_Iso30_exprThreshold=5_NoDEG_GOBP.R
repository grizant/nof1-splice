## Score TCGA UCEC patients using 2018 GOBP
## Remove DEG from edgeR analysis and low expressing genes 
## AG Schissler
## Created 9 Aug 2018

##############################################################################
#### 1. Setup environment

## DO NOT RUN (here for reproducibility)
## create GO-BP from Jianrong's 2018 gene set definitions
## gobp <- read.delim("~/Dropbox/Lab-Tools-Files/GeneOntology/2018-02-07/Human/gene2go_Human_BP_filter.txt")
## head(gobp)
## length(unique(gobp$goid))
## sizes <- unlist(lapply(split(gobp$symbol, gobp$goid), length))
## str(sizes); summary(sizes)
## genes <- unique(as.character(gobp$symbol))
## str(genes) ## 15290 in question

## filter UCEC iso data to just genes in pathways
## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/ucec_iso_paired.RData")
## ucec_iso_gobp_data <- ucec_iso[which(ucec_iso$geneSymbol %in% genes), ]
## str(ucec_iso_gobp_data)
## tmp_genes <- unique(as.character(gobp$symbol))
## saveRDS(ucec_iso_gobp_data, file = "~/Dropbox/Splice-n-of-1-pathways/Data/ucec_iso_gobp_data.rds")

## load TCGA UCEC TPM isoform data
ucec_iso_gobp_data <- readRDS(file = "~/Dropbox/Splice-n-of-1-pathways/Data/ucec_iso_gobp_data.rds") 
## DEG isoform data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/ucec_deg_list_all.RData") 

## source functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R")

##############################################################################
#### 2. Restructure iso data into a patient-wise list for parallel processing

## Retrieve patient IDs
pat_col <- grep("TCGA", x = names(ucec_iso_gobp_data))
patients_chr <- unique(substring(names(ucec_iso_gobp_data[pat_col]), 1, 12))

## create a empty list
iso_gobp_list <- vector(mode = "list", length =  length(patients_chr))
names(iso_gobp_list) <- patients_chr

## tmp_pat <- patients_chr[1]
for (tmp_pat in patients_chr) {
    ## retrieve gene symbols and the paired transcriptomes
    iso_gobp_list[[tmp_pat]] <- (data.frame(geneSymbol = ucec_iso_gobp_data$geneSymbol,
                                            ucec_iso_gobp_data[, grep(tmp_pat, names(ucec_iso_gobp_data))]))
    iso_gobp_list[[tmp_pat]][, "geneSymbol"] <- as.character(ucec_iso_gobp_data$geneSymbol)
}

##############################################################################
#### 2b. Restructure DEG data into a patient-wise list for parallel processing

## DEG_dat <- ucec_deg_list[[1]]

get_DEGs <- function(DEG_dat, threshold = 0.05) {
    DEG_dat$id[DEG_dat$BY <= threshold]
}

## 20% to follow the rest of the paper
## this is lenient and will remove more genes
DEG_list <- lapply(ucec_deg_list_all, get_DEGs, threshold=0.2)

##############################################################################
#### 3. Setup parallel processing

## load parallelization library
library(parallel)

## start a cluster
num_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(num_cores, type = "FORK")

## export data and objects

## export custom functions
parallel::clusterEvalQ(cl, expr = source("~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R"))

## export local FDR of Efron
parallel::clusterEvalQ(cl, expr = library(locfdr))

## run on a small subset
## set.seed(44444)
## tmp_index <- sample(1:length(patients_chr), size = num_cores)
## ## 
## tmp_list <- iso_gobp_list[tmp_index]
## tmp_deg_list <- DEG_list[tmp_index]
## 
## ## use mapply since we have to change inputs per patient
## avg_scores <- parallel::mcmapply(function(X, Y) {
##     transform_iso_pathway(iso_data=X, DEGs=Y, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GOBP/gobp_tb.txt",
## desc_file = "~/Dropbox/Lab-Tools/GeneSets/GOBP/gobp.description_tb.txt", pathway_method = "avg", gene_method = "hellinger")
## }, X = tmp_list, Y = DEG_list, SIMPLIFY = F)
## str(avg_scores)
 
## system.time(avg_scores <- parallel::parLapply(cl = cl, tmp_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GOBP/gobp_tb.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/GOBP/gobp.description_tb.txt", pathway_method = "avg", gene_method = "hellinger"))
## 
## scores_list <- avg_scores

##############################################################################
#### 4. Score all patients

system.time(scores_list <- parallel::mcmapply(function(X, Y) {
    transform_iso_pathway(iso_data=X, DEGs=Y, annot_file = "~/Dropbox/Lab-Tools-Files/GeneOntology/2018-02-07/Human/gene2go_Human_BP_filter.txt",
desc_file = "~/Dropbox/Lab-Tools-Files/GeneOntology/2018-02-07/gene_ontology02072018.txt", pathway_method = "EE", gene_method = "hellinger", genes_range = c(15,500), expr_threshold = 5, go_desc_file = T, path_id_name = "goid", database = "gobp")
}, X = iso_gobp_list, Y = DEG_list, SIMPLIFY = F))

## str(scores_list)
## save the object
save(scores_list, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_UCEC_hel_EE_Iso30_exprThreshold=5_NoDEG_GOBP_9aug2018.RData")

## close cluster
parallel::stopCluster(cl = cl)

##############################################################################
#### 5. Explore quickly

(num_hits <- unlist(lapply(scores_list, function(tmp_data){
    sum(tmp_data$diff_splice_call, na.rm = T)
})))

summary(num_hits)

lapply(scores_list, function(tmp_data){
    summary(tmp_data$fdr_value)
    head(tmp_data, 10)
})
