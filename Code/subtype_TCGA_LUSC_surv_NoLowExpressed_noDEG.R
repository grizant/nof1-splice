## LUSC subtyping via Splice-N-of-1-pathways
## AG Schissler
## Created 27 Jul 2017

############################################################
## i. Load clinical data and results

load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_paired_clinical.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expr_threshold=5_NoDEG_KEGG_25july2018.RData")
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
threshold <- 0.20
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
top_hits <- names(transformed_pvalue)[transformed_pvalue < tmp_locfdr$z.2[1]]
tmp_data[top_hits,] ## top hits are all interesting including staph infection!

######################################################
## explore top hits


### hit 1
i <- 1
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p1 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p1 <- p1 + geom_boxplot())

### hit 2
i <- 2
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p2 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p2 <- p2 + geom_boxplot())

### hit 3
i <- 3
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p3 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p3 <- p3 + geom_boxplot())

### hit 4
i <- 4
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p4 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p4 <- p4 + geom_boxplot())
