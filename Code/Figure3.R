## Create Figure 3
## AG Schissler
## Created 27 Jul 2017

############################################################
## i. Set up environment

source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R")
source("~/Dropbox/Splice-n-of-1-pathways/Code/target_functions.R")
library(ggplot2)
library(cowplot)

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")

############################################################
## 1. Organize target pathway metric in one data frame

## load and aggregate manually
all_data <- NULL

## 1. BLCA
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(cancer_pathways <- get_cancer_pathways(scores_list))
(cancer_ids <- names(cancer_pathways))
(target <- names(cancer_pathways)[grep("Bladder", cancer_pathways)])
(target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T))
target_data$TCGA <- "BLCA"
## add cluster to show patient heterogeneity
target_data$cluster <- cluster_pat(target_data$OR)
(all_data <- rbind(all_data, target_data))

## 2. THCA
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_THCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Thyroid", cancer_pathways)])
(target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T))
target_data$TCGA <- "THCA"
## add cluster to show patient heterogeneity
target_data$cluster <- cluster_pat(target_data$OR)
(all_data <- rbind(all_data, target_data))

## 3. UCEC
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_UCEC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Endo", cancer_pathways)])
(target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T))
target_data$TCGA <- "UCEC"
## add cluster to show patient heterogeneity
target_data$cluster <- cluster_pat(target_data$OR)
(all_data <- rbind(all_data, target_data))

## 4. PRAD
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_PRAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Prostate", cancer_pathways)])
(target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T))
target_data$TCGA <- "PRAD"
## add cluster to show patient heterogeneity
target_data$cluster <- cluster_pat(target_data$OR)
(all_data <- rbind(all_data, target_data))

## 5. LUSC
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
(target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T))
target_data$TCGA <- "LUSC"
## add cluster to show patient heterogeneity
sort(target_data$OR)
(target_data$cluster <- cluster_pat(target_data$OR))
sort(target_data$OR[target_data$cluster == 2])
(all_data <- rbind(all_data, target_data))

## 6. LUAD
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
(target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T))
target_data$TCGA <- "LUAD"
## add cluster to show patient heterogeneity
target_data$cluster <- cluster_pat(target_data$OR)
(all_data <- rbind(all_data, target_data))

############################################################
## 2. create plot

str(all_data)
nrow(all_data)

plot_data <- all_data
### factor management
plot_data$TCGA <- factor(plot_data$TCGA)
## levels(plot_data$TCGA)
## plot_data$TCGA <- factor(plot_data$TCGA, levels = c("LUSC", "LUAD", "PRAD", "THCA", "UCEC", "BLCA"))
## plot_data$TCGA <- factor(plot_data$TCGA, labels = c("BLCA (19)", "LUAD (58)", "LUSC (51)", "PRAD (52)", "THCA (59)", "UCEC (7)"))

table(plot_data$TCGA)
sort(plot_data[plot_data$TCGA == "LUSC", "OR"])
plot_data$cluster <- factor(plot_data$cluster)
library(dplyr)
grouped <- group_by(plot_data, TCGA, cluster)
summarise(grouped, mean=mean(OR), sd=sd(OR))

## create interaction
## plot_data$inter <- factor(paste(as.character(plot_data$TCGA), as.character(plot_data$cluster), sep="-"), labels = as.character(rep(1:2, times = length(unique(plot_data$TCGA)))))
plot_data$inter <- paste(as.character(plot_data$TCGA), as.character(plot_data$cluster), sep="-")

## levels(plot_data$inter) <- levels(plot_data$inter)[c(1,7,2,8,3,9,4,10,5,11,6,12)]

(fig3 <- ggplot(data = plot_data, aes(x = inter, y = OR, color = TCGA)) +
     geom_boxplot() +
     theme(legend.position = "right") +
     background_grid(major = "xy", minor = "none") +
     ylim(c(0, 4)) +
     geom_jitter(width = 0.1, alpha = 0.5) +
     geom_hline(yintercept = 1) +
     xlab( "Cluster via KEGG target pathway" ) +
     ylab( "Odds ratio of enriched alt splice") +
     scale_color_discrete("TCGA", 
                         labels = c("BLCA (19)", "LUAD (58)", "LUSC (51)", "PRAD (52)", "THCA (59)", "UCEC (7)"))
    )

save_plot("Figure3.pdf", fig3, base_height = 8, base_aspect_ratio = 1.3)

