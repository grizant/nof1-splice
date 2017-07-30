## Create Table 4
## AG Schissler
## Created 27 Jul 2017

############################################################
## i. Load clinical data and results

load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_paired_clinical.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
patients <- names(scores_list)
clin_data <- clin_data[rownames(clin_data) %in% patients,]

## load survival analysis functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R")

## explore clinical
table(clin_data$vitalstatus)

############################################################
## 1. Aggregate pathway scores

effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)

## filter to only those pathways found dysregulated in at least one patient
fil_effect_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "pathway_score")
fil_fdr_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "fdr_value")
str(fil_effect_mat) ## 100 pathways now
str(fil_fdr_mat)

############################################################
## 2. Select pathways that produce distinct survival curves

## two clusters (not sure how to justify this other than simplicity)
clust_pvalue <- get_pathway_pvalue(value_mat = fil_effect_mat, clin_data = clin_data, type = "pam", num_clusters = 2)
hist(clust_pvalue, breaks = 15)
clust_pvalue <- sort(clust_pvalue)
head(clust_pvalue)

## explore hits
threshold <- 0.05
any(clust_pvalue < threshold)
(hits <- names(clust_pvalue)[clust_pvalue < threshold])
sum(clust_pvalue < threshold)/length(clust_pvalue)

## retrive the hit pathway description 
tmp_data <- scores_list[[1]]
tmp_data[hits,]

#### try fdr adjustments
## standard adjustment
## p.adjust(p = clust_pvalue, "fdr")

## locFDR works and gets around the permutation approach!
transformed_pvalue <- qnorm(clust_pvalue)
eps <- 0.001
transformed_pvalue[is.infinite(transformed_pvalue)] <- max(transformed_pvalue[!is.infinite(transformed_pvalue)], na.rm = T) + eps

tmp_locfdr <- locfdr::locfdr(zz = transformed_pvalue, bre = ceiling(length(clust_pvalue)/8), df = 4, pct = 0, pct0 = 1/64, nulltype = 1, plot = 1)
tmp_locfdr

## 4 hits
head(tmp_locfdr$fdr)
top_hits <- names(transformed_pvalue)[transformed_pvalue < tmp_locfdr$z.2[1]]
tmp_data[top_hits,] ## top hits are all interesting including staph infection!

###### create Table 4

## number of genes
tmp_data <- scores_list[[1]]
tmp_data[top_hits,]

## gene overlap
kegg <- nof1::read_gene_set("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")
kegg_list <- split(kegg$symbol, kegg$path_id)
top_genes <- kegg_list[top_hits]
tmp_genes <- top_genes[[2]]
target_genes <- top_genes[[1]]
lapply(top_genes, function(tmp_genes) {
    length(intersect(x = tmp_genes, y = target_genes))/length(tmp_genes)
})
    
## percent calls
top_hit <- top_hits[1]

for (top_hit in top_hits) {
    print(sum(unlist(lapply(scores_list, function(tmp_data) {
        tmp_data[top_hit, "diff_splice_call"]
    })))/length(scores_list)*100)
}
