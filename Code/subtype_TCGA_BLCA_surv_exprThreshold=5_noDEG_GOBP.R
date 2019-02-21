## BLCA subtyping via Splice-N-of-1-pathways on GOBP pathways
## AG Schissler
## Created 10 Aug 2018

############################################################
## i. Load clinical data and results

load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_paired_clinical.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_exprThreshold=5_NoDEG_GOBP_9aug2018.RData")
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
str(fil_effect_mat) ## 1180 pathways now
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
## blarg. looks uniform... with almost the exact expected proportion

## retrive the hit pathway description 
tmp_data <- scores_list[[1]]
tmp_data[hits,]

#### try fdr adjustments
## standard adjustment
clust_adj_pvalue <- p.adjust(p = clust_pvalue, "fdr")
summary(clust_adj_pvalue)

## ## locFDR works and gets around the permutation approach!
## transformed_pvalue <- qnorm(clust_pvalue)
## eps <- 0.001
## transformed_pvalue[is.infinite(transformed_pvalue)] <- max(transformed_pvalue[!is.infinite(transformed_pvalue)], na.rm = T) + eps
## 
## tmp_locfdr <- locfdr::locfdr(zz = transformed_pvalue, plot = 1)
## tmp_locfdr
## 
## 
## top_hits <- names(transformed_pvalue)[transformed_pvalue < tmp_locfdr$z.2[1]]
top_hits <- names(which(clust_adj_pvalue < threshold))
tmp_data[names(which(clust_adj_pvalue < threshold)),]

## 2 hits

######################################################
## explore top hits


### hit 1 (nope only 1 is the other cluster)
i <- 1
(top_or <- fil_effect_mat[,top_hits[i]])
sort((top_fdr <- fdr_mat[,top_hits[i]]))
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p1 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p1 <- p1 + geom_boxplot())

### hit 2 (same)
i <- 2
(top_or <- fil_effect_mat[,top_hits[i]])
sort((top_fdr <- fdr_mat[,top_hits[i]]))
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
