## Score TCGA THCA patients via Average pathway Hellinger distances
## and Empirical Enrichment
## Remove DEG from edgeR analysis
## AG Schissler
## Created 25 Jul 2018

##############################################################################
#### 1. Setup environment

## load TCGA THCA TPM isoform data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/thca_iso_kegg_data.RData") ## NEW ISOFORM 25 Jul 2017
## DEG isoform data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/thca_deg_list_all.RData") 

## source functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R")

##############################################################################
#### 2. Restructure iso data into a patient-wise list for parallel processing

## Retrieve patient IDs
pat_col <- grep("TCGA", x = names(thca_iso_kegg_data))
patients_chr <- unique(substring(names(thca_iso_kegg_data[pat_col]), 1, 12))

## create a empty list
iso_kegg_list <- vector(mode = "list", length =  length(patients_chr))
names(iso_kegg_list) <- patients_chr

## tmp_pat <- patients_chr[1]
for (tmp_pat in patients_chr) {
    ## retrieve gene symbols and the paired transcriptomes
    iso_kegg_list[[tmp_pat]] <- (data.frame(geneSymbol = thca_iso_kegg_data$geneSymbol,
                                            thca_iso_kegg_data[, grep(tmp_pat, names(thca_iso_kegg_data))]))
    iso_kegg_list[[tmp_pat]][, "geneSymbol"] <- as.character(thca_iso_kegg_data$geneSymbol)
}

##############################################################################
#### 2b. Restructure DEG data into a patient-wise list for parallel processing

## DEG_dat <- thca_deg_list[[1]]

get_DEGs <- function(DEG_dat, threshold = 0.05) {
    DEG_dat$id[DEG_dat$BY <= threshold]
}

## 20% to follow the rest of the paper
## this is lenient and will remove more genes
DEG_list <- lapply(thca_deg_list_all, get_DEGs, threshold=0.2)

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
set.seed(44444)
tmp_index <- sample(1:length(patients_chr), size = num_cores)
## 
tmp_list <- iso_kegg_list[tmp_index]
tmp_deg_list <- DEG_list[tmp_index]

## use mapply since we have to change inputs per patient
avg_scores <- parallel::mcmapply(function(X, Y) {
    transform_iso_pathway(iso_data=X, DEGs=Y, annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt",
desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", pathway_method = "avg", gene_method = "hellinger")
}, X = tmp_list, Y = DEG_list, SIMPLIFY = F)
str(avg_scores)

## system.time(avg_scores <- parallel::parLapply(cl = cl, tmp_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", pathway_method = "avg", gene_method = "hellinger"))
## 
## scores_list <- avg_scores

##############################################################################
#### 4. Score all patients

## lower the minimum number of genes for the w/o DEG analysis
## (important pathways are being excluded)
## now score by Empirical Enrichment (~9.4 sec)
system.time(scores_list <- parallel::mcmapply(function(X, Y) {
    transform_iso_pathway(iso_data=X, DEGs=Y, annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt",
desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", pathway_method = "EE", gene_method = "hellinger", genes_range = c(5,500))
}, X = iso_kegg_list, Y = DEG_list, SIMPLIFY = F))

## str(scores_list)
## save the object
save(scores_list, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_THCA_hel_EE_Iso30_expressed_NoDEG_KEGG_25july2018.RData")

## close cluster
parallel::stopCluster(cl = cl)

##############################################################################
#### 5. Explore quickly

(num_hits <- unlist(lapply(scores_list, function(tmp_data){
    sum(tmp_data$diff_splice_call, na.rm = T)
})))

summary(num_hits)

## find number of cancer hits
## tmp_data <- scores_list[[1]]
## tmp_id <- "hsa05200"

(cancer_logic <- unlist(lapply(scores_list, function(tmp_data){
    ## tmp_data[tmp_id, ]
    any(grepl("cancer", tmp_data[tmp_data$diff_splice_call == 1, "pathway_desc"]))
})))

## tmp_data[grep("cancer", tmp_data$pathway_desc),]
any(cancer_logic)
sum(cancer_logic)/length(cancer_logic)

summary(num_hits)
## qplot(num_hits)

## explore some results
set.seed(44)
tmp_data <- scores_list[[sample(1:length(scores_list),1)]]
## str(tmp_data)
## head(tmp_data, 20)
tmp_data[grep("cancer", tmp_data$pathway_desc),]
tmp_data[grep("Bladder cancer", tmp_data$pathway_desc),]

## find rank of Bladder cancer pathway, hsa05219
tmp_id <- "hsa05219"

unlist(lapply(scores_list, function(tmp_data) {
    which(rownames(tmp_data) == tmp_id)
}))

unlist(lapply(scores_list, function(tmp_data) {
    tmp_data[tmp_id, "pathway_score"]
}))

tmp_pat <- "TCGA-BT-A20U"
nrow(scores_list[[tmp_pat]])
scores_list[[tmp_pat]][165,]

##############################################################################
#### 6. Systematically explore

load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_THCA_hel_EE_Iso30_expressed_NoDEG_KEGG_25july2018.RData")

## 6.1. Capture rate of target pathway while varying FDR
## find rank of Bladder cancer pathway, hsa05219
target_id <- "hsa05219"

## tmp_data <- scores_list[[1]]

lapply(scores_list, function(tmp_data){
    summary(tmp_data$fdr_value)
    head(tmp_data, 10)
})
