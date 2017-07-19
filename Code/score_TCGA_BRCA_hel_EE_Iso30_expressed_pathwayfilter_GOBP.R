## Score TCGA BRCA patients via Empirical Enrichment of Hellinger distances
## AG Schissler
## Created 18 Jul 2017

##############################################################################
#### 1. Setup environment

## load TCGA BRCA TPM isoform data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_tpm_GO-BP.RData")

## source functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R")

##############################################################################
#### 2. Restructure iso data into a patient-wise list for parallel processing

## Retrieve patient IDs
patients_chr <- unique(substring(names(iso_go_data[,-(1:4)]), 1, 12))

## create a empty list
iso_go_list <- vector(mode = "list", length =  length(patients_chr))
names(iso_go_list) <- patients_chr

## tmp_pat <- patients_chr[1]
for (tmp_pat in patients_chr) {
    ## retrieve gene symbols and the paired transcriptomes
    iso_go_list[[tmp_pat]] <- (data.frame(Gene_symbol = iso_go_data$Gene_symbol,
                                          iso_go_data[, grep(tmp_pat, names(iso_go_data))]))
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
## set.seed(444)
## tmp_index <- sample(1:length(patients_chr), size = num_cores)
## 
## tmp_list <- iso_go_list[tmp_index]
## 
## system.time(avg_scores <- parallel::parLapply(cl = cl, tmp_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_description.txt", pathway_method = "EE", gene_method = "hellinger"))
## 55 seconds
## scores_list <- avg_scores

##############################################################################
#### 4. Score all patients

system.time(scores_list <- parallel::parLapply(cl = cl, iso_go_list, transform_iso_pathway, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_description.txt", pathway_method = "EE", gene_method = "hellinger")) ## 390 seconds

## save the object
save(scores_list, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_EE_Iso30_expressiod_pathwayfilter_13july2017.RData")

## close cluster
parallel::stopCluster(cl = cl)

##############################################################################
#### 5. Explore

load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_avg_Iso30_expressiod_pathwayfilter_13july2017.RData")

(num_hits <- unlist(lapply(scores_list, function(tmp_data){
    sum(tmp_data$diff_splice_call, na.rm = T)
})))

summary(num_hits)
qplot(num_hits)

## isoform distance relationship

tmp_iso <- iso_go_list[[1]]

## find average number of isoforms 
str(tmp_iso)
