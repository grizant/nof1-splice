## Validate TCGA BLCA patients via Average pathway Hellinger distances
## and Empirical Enrichment
## applying filters to reduce noise
## AG Schissler
## Created 25 Jul 2017

##############################################################################
#### 1. Load and systematically explore

load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")

## 6.1. Capture rate of target pathway while varying FDR
## find rank of Bladder cancer pathway, hsa05219
target_id <- "hsa05219"


## tmp_data <- scores_list[[1]]

lapply(scores_list, function(tmp_data){
    summary(tmp_data$fdr_value)
    head(tmp_data$fdr_value, 10)
})
