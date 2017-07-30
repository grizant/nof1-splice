## Create Figures 4 and 5
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

kegg <- nof1::read_gene_set(file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")
kegg_desc <- read.delim("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", stringsAsFactors = F)
rownames(kegg_desc) <- as.character(kegg_desc$path_id)

## explore clinical
table(clin_data$vitalstatus)

############################################################
## 1. Aggregate pathway scores

effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)

## filter to only those pathways found dysregulated in at least one patient
fil_effect_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "pathway_score")
fil_fdr_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "fdr_value")
str(fil_effect_mat) ## 101 pathways now
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
kegg_desc[hits,]

#### try fdr adjustments
## standard adjustment
## p.adjust(p = clust_pvalue, "fdr")

## locFDR works and gets around the permutation approach!
hist(clust_pvalue, breaks = 30)
transformed_pvalue <- qnorm(clust_pvalue)
eps <- 0.0001
transformed_pvalue[is.infinite(transformed_pvalue)] <- max(transformed_pvalue[!is.infinite(transformed_pvalue)], na.rm = T) + eps
hist(transformed_pvalue, breaks = 30)

tmp_locfdr <- locfdr::locfdr(zz = clust_pvalue, bre = ceiling(length(clust_pvalue)/8), df = 4, pct = 0, pct0 = 1/4, nulltype = 1, plot = 1)

tmp_locfdr <- locfdr::locfdr(zz = transformed_pvalue, bre = ceiling(length(clust_pvalue)/8), df = 4, pct = 0, pct0 = 1/4, nulltype = 1, plot = 1)

sum(transformed_pvalue < tmp_locfdr$z.2[1])

## 4 hits
threshold <- 0.5
top_hits <- names(transformed_pvalue)[transformed_pvalue < tmp_locfdr$z.2[1]]
top_hits <- names(transformed_pvalue)[tmp_locfdr$fdr < threshold]
tmp_data[top_hits,] ## top hits are all interesting including staph infection!
tmp_locfdr$fdr[1:length(top_hits)]

hits <- names(transformed_pvalue)[1:4]
kegg_desc[hits,]

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
(top_desc <- tmp_data[hits[i], "pathway_desc"])
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
(top_desc <- tmp_data[hits[i], "pathway_desc"])
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
(top_desc <- tmp_data[hits[i], "pathway_desc"])
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
top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T)
target_clust <- four_clust
(top_desc <- tmp_data[hits[i], "pathway_desc"])
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

## see agreement in clustering

for (tmp_clust in list(one_clust, two_clust, three_clust)) {
    print(sum(target_clust == tmp_clust)/length(target_clust)*100)
}

##################
### combine together

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")

## Figure 4
surv_plots <- plot_grid(s1, s2, s3, s4,
                         labels=c("A", "B", "C", "D"), ncol = 2, align = "h")
save_plot("Figure4.pdf", surv_plots,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.1,
          )

## Figure 5
or_plots <- plot_grid(p1, p2, p3, p4,
                         labels=c("A", "B", "C", "D"), ncol = 2, align = "h")
save_plot("Figure5.pdf", or_plots,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.1,
          base_width = 5
          )
