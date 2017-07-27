## Validate TCGA COAD patients via Average pathway Hellinger distances
## and Empirical Enrichment
## applying filters to reduce noise
## AG Schissler
## Created 25 Jul 2017

##############################################################################
#### 1. Load and systematically explore

load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_COAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")

source("~/Dropbox/Splice-n-of-1-pathways/Code/target_functions.R")

(cancer_pathways <- get_cancer_pathways(scores_list))
(cancer_ids <- names(cancer_pathways))
(target <- names(cancer_pathways)[grep("Colorectal", cancer_pathways)])


target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
summarize_target_data(target_data)
