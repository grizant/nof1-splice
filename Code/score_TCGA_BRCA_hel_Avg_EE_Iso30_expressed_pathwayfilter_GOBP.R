## Score TCGA BRCA patients via Average pathway Hellinger distances
## and Empirical Enrichment
## applying filters to reduce noise
## AG Schissler
## Created 25 Jul 2017

##############################################################################
#### 1. Setup environment

## load TCGA BRCA TPM isoform data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/brca_iso_gobp_data.RData") ## NEW ISOFORM 25 Jul 2017

## source functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R")

##############################################################################
#### 2. Restructure iso data into a patient-wise list for parallel processing

## Retrieve patient IDs
patients_chr <- unique(substring(names(brca_iso_gobp_data[,-(1:4)]), 1, 12))

## create a empty list
iso_gobp_list <- vector(mode = "list", length =  length(patients_chr))
names(iso_gobp_list) <- patients_chr

## tmp_pat <- patients_chr[1]
for (tmp_pat in patients_chr) {
    ## retrieve gene symbols and the paired transcriptomes
    iso_gobp_list[[tmp_pat]] <- (data.frame(geneSymbol = brca_iso_gobp_data$geneSymbol,
                                            brca_iso_gobp_data[, grep(tmp_pat, names(brca_iso_gobp_data))]))
    iso_gobp_list[[tmp_pat]][, "geneSymbol"] <- as.character(brca_iso_gobp_data$geneSymbol)
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
tmp_list <- iso_gobp_list[tmp_index]
system.time(avg_scores <- parallel::parLapply(cl = cl, tmp_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_description.txt", pathway_method = "EE", gene_method = "hellinger"))
## 25 seconds
## scores_list <- avg_scores

##############################################################################
#### 4. Score all patients

system.time(scores_list <- parallel::parLapply(cl = cl, iso_gobp_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_description.txt", pathway_method = "avg", gene_method = "hellinger")) ## 196.8 seconds

## save the object
save(scores_list, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_avg_Iso30_expressiod_pathwayfilter_GOBP_25july2017.RData")

## now score by Empirical Enrichment
system.time(scores_list <- parallel::parLapply(cl = cl, iso_gobp_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_description.txt", pathway_method = "EE", gene_method = "hellinger")) ## 426 seconds

## save the object
save(scores_list, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_EE_Iso30_expressiod_pathwayfilter_GOBP_25july2017.RData")

## close cluster
parallel::stopCluster(cl = cl)

##############################################################################
#### 5. Explore

(num_hits <- unlist(lapply(scores_list, function(tmp_data){
    sum(tmp_data$diff_splice_call, na.rm = T)
})))

summary(num_hits)

## explore some results
set.seed(44)
tmp_data <- scores_list[[sample(1:length(scores_list),1)]]
str(tmp_data)
head(tmp_data, 20)


