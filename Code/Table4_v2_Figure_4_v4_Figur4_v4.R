## Create Table 4, locFDR was unstable try new attempt
## AG Schissler
## Created 30 Jul 2017

############################################################
## i. Load clinical data and results

load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_paired_clinical.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
patients <- names(scores_list)
clin_data <- clin_data[rownames(clin_data) %in% patients,]

## load survival analysis functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R")

kegg <- nof1::read_gene_set(file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")
kegg_desc <- read.delim("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", stringsAsFactors = F)
rownames(kegg_desc) <- as.character(kegg_desc$path_id)

## explore clinical
table(clin_data$vitalstatus)

############################################################
## 1. Aggregate pathway scores

effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
rm_effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = T)
str(rm_effect_mat)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)

## filter to only those pathways found dysregulated in at least one patient
fil_effect_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "pathway_score")
fil_fdr_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "fdr_value")
str(fil_effect_mat) ## 101 pathways now
str(fil_fdr_mat)

############################################################
## 2. Select pathways that produce distinct survival curves

## get empirical distribution by sampling
set.seed(444)
reps <- 2000
## tmp_pat <- fil_effect_mat[1,]
rand_effect_mat <- NULL
for (i in 1:reps) {
    rand_effect_mat <- cbind(rand_effect_mat, apply(fil_effect_mat, 1, function(tmp_pat) {
        unname(sample(tmp_pat, 1))
    }))
}

## two clusters (not sure how to justify this other than simplicity)
## clust_pvalue <- get_pathway_pvalue(value_mat = rm_effect_mat, clin_data = clin_data, type = "pam", num_clusters = 2)
clust_pvalue <- get_pathway_pvalue(value_mat = fil_effect_mat, clin_data = clin_data, type = "pam", num_clusters = 2)
clust_pvalue <- sort(clust_pvalue)

rand_pvalue <- get_pathway_pvalue(value_mat = rand_effect_mat, clin_data = clin_data, type = "pam", num_clusters = 2)

emp_pvalue <- sapply(clust_pvalue, function(tmp_pvalue) {
    tmp_p <- sum(rand_pvalue < tmp_pvalue)/reps
    if (tmp_p == 0) tmp_p <- (1/(reps + 1))
    tmp_p
})

emp_pvalue <- sort(emp_pvalue)
head(emp_pvalue, 6)
emp_fdr <- p.adjust(emp_pvalue, method = "fdr")
head(emp_fdr, 6)

## explore hits
threshold <- 0.20
any(emp_fdr < threshold)
(hits <- names(clust_pvalue)[which(emp_fdr <= threshold)])

## retrive the hit pathway description 
kegg_desc[hits,]
head(clust_pvalue, length(hits)+1)
head(emp_pvalue, length(hits)+1)
head(emp_fdr, length(hits)+1)

######################################################
## explore top hits
library(ggplot2)
library(cowplot)

### hit 1
i <- 1
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
one_clust <- top_clust
(top_desc <- kegg_desc[hits[i], "description"])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s1 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))
## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Better", "Worse")))
(p1 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio)) +
     geom_boxplot() + 
     ## theme_bw() +
     theme(legend.position = "none") +
     background_grid(major = "xy", minor = "none") +
     ylim(c(0, 3.25)) +
     geom_jitter(width = 0.1, alpha = 0.5) +
     geom_hline(yintercept = 1) +
     xlab( paste("Cluster via '", top_desc, "'", sep="") ) )
    
### hit 2
i <- 2
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
two_clust <- top_clust
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(top_desc <- kegg_desc[hits[i], "description"])
(s2 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))
## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Better", "Worse")))
(p2 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio)) +
     geom_boxplot() + 
     ## theme_bw() +
     theme(legend.position = "none") +
     background_grid(major = "xy", minor = "none") +
     ylim(c(0, 3.25)) +
     geom_jitter(width = 0.1, alpha = 0.5) +
     geom_hline(yintercept = 1) +
     xlab( paste("Cluster via '", top_desc, "'", sep="") ) )

### hit 3
i <- 3
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
three_clust <- top_clust
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(top_desc <- kegg_desc[hits[i], "description"])
top_desc <- gsub("Staphylococcus", "Staph", top_desc)
(s3 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))
## cluster 2 is better!!
top_clust <- ifelse(top_clust == 1, 2, 1)
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Better", "Worse")))
(p3 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio)) +
     geom_boxplot() + 
     ## theme_bw() +
     theme(legend.position = "none") +
     background_grid(major = "xy", minor = "none") +
     ylim(c(0, 3.25)) +
     geom_jitter(width = 0.1, alpha = 0.5) +
     geom_hline(yintercept = 1) +
     xlab( paste("Cluster via '", top_desc, "'", sep="") ) )

### hit 4
i <- 4
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
four_clust <- top_clust
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
target_clust <- four_clust
(top_desc <- kegg_desc[hits[i], "description"])
## top_desc <- "Staph aureus infection"
(s4 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))
## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Better", "Worse")))
(p4 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio)) +
     geom_boxplot() + 
     ## theme_bw() +
     theme(legend.position = "none") +
     background_grid(major = "xy", minor = "none") +
     ylim(c(0, 3.25)) +
     geom_jitter(width = 0.1, alpha = 0.5) +
     geom_hline(yintercept = 1) +
     xlab( paste("Cluster via '", top_desc, "'", sep="") ) )

### hit 5
i <- 5
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
five_clust <- top_clust
top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T)
top_surv
(top_desc <- kegg_desc[hits[i], "description"])
## top_desc <- "Staph aureus infection"
(s5 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))
## cluster 2 is better
top_clust <- ifelse(top_clust == 1, 2, 1)
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Better", "Worse")))
(p5 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio)) +
     geom_boxplot() + 
     ## theme_bw() +
     theme(legend.position = "none") +
     background_grid(major = "xy", minor = "none") +
     ylim(c(0, 3.25)) +
     geom_jitter(width = 0.1, alpha = 0.5) +
     geom_hline(yintercept = 1) +
     xlab( paste("Cluster via '", top_desc, "'", sep="") ) )


## see agreement in clustering

## sum(three_clust == five_clust)/length(three_clust)

for (tmp_clust in list(one_clust, two_clust, three_clust, four_clust, five_clust)) {
    print(sum(target_clust == tmp_clust)/length(target_clust)*100)
}

##################
### combine together

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")

## Figure 4
(surv_plots <- plot_grid(s1, s2, s3, s4,
                         labels=c("A", "B", "C", "D"), ncol = 2, align = "h"))
save_plot("Figure4_v4.pdf", surv_plots,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.1,
          )

## Figure 5
(or_plots <- plot_grid(p1, p2, p3, p4,
                         labels=c("A", "B", "C", "D"), ncol = 2, align = "h"))
save_plot("Figure5_v4.pdf", or_plots,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.1,
          base_width = 5
          )

###### create more summaries forTable 4

## number of genes

## gene overlap
kegg <- nof1::read_gene_set("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")
kegg_list <- split(kegg$symbol, kegg$path_id)
top_genes <- kegg_list[hits]
tmp_genes <- top_genes[[1]]
length(top_genes[[5]])
target_genes <- top_genes[["hsa05223"]]
lapply(top_genes, function(tmp_genes) {
    length(intersect(x = tmp_genes, y = target_genes))/length(tmp_genes)
})
    
## percent calls


for (top_hit in hits) {
    print(sum(unlist(lapply(scores_list, function(tmp_data) {
        tmp_data[top_hit, "diff_splice_call"]
    })))/length(scores_list)*100)
}
