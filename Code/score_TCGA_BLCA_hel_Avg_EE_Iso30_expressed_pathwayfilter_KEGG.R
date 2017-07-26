## Score TCGA BLCA patients via Average pathway Hellinger distances
## and Empirical Enrichment
## applying filters to reduce noise
## AG Schissler
## Created 25 Jul 2017

##############################################################################
#### 1. Setup environment

## load TCGA BLCA TPM isoform data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_kegg_data.RData") ## NEW ISOFORM 25 Jul 2017

## source functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R")

##############################################################################
#### 2. Restructure iso data into a patient-wise list for parallel processing

## Retrieve patient IDs
patients_chr <- unique(substring(names(blca_iso_kegg_data[,-(1:4)]), 1, 12))

## create a empty list
iso_kegg_list <- vector(mode = "list", length =  length(patients_chr))
names(iso_kegg_list) <- patients_chr

## tmp_pat <- patients_chr[1]
for (tmp_pat in patients_chr) {
    ## retrieve gene symbols and the paired transcriptomes
    iso_kegg_list[[tmp_pat]] <- (data.frame(geneSymbol = blca_iso_kegg_data$geneSymbol,
                                            blca_iso_kegg_data[, grep(tmp_pat, names(blca_iso_kegg_data))]))
    iso_kegg_list[[tmp_pat]][, "geneSymbol"] <- as.character(blca_iso_kegg_data$geneSymbol)
}

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
system.time(avg_scores <- parallel::parLapply(cl = cl, tmp_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", pathway_method = "EE", gene_method = "hellinger"))
## 
## scores_list <- avg_scores

##############################################################################
#### 4. Score all patients

system.time(scores_list <- parallel::parLapply(cl = cl, iso_kegg_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", pathway_method = "avg", gene_method = "hellinger")) ## 4 seconds

## save the object
save(scores_list, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_avg_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")

## now score by Empirical Enrichment
system.time(scores_list <- parallel::parLapply(cl = cl, iso_kegg_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", pathway_method = "EE", gene_method = "hellinger")) ## 5 seconds

## save the object
save(scores_list, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")

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
qplot(num_hits)

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
    tmp_data[tmp_id, "pathway_score"
}))

tmp_pat <- "TCGA-BT-A20U"
nrow(scores_list[[tmp_pat]])
scores_list[[tmp_pat]][165,]

##############################################################################
#### 6. Systematically explore

load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")

## 6.1. Capture rate of target pathway while varying FDR
## find rank of Bladder cancer pathway, hsa05219
target_id <- "hsa05219"


## tmp_data <- scores_list[[1]]

lapply(scores_list, function(tmp_data){
    summary(tmp_data$fdr_value)
    head(tmp_data$fdr_value, 10)
})
