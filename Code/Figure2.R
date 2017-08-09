## Create Figure 2
## AG Schissler
## Created 30 Jul 2017

############################################################
## i. Set up environment

library(ggplot2)
library(cowplot)
library(png)
library(grid)

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure1_dist_data.RData")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure2_emp_pvalue.RData")
load(file = "~/Dropbox/splice-n-of-1-pathways/Data/Figure1_odds_ratio.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/lusc_iso_kegg_data.RData") ## NEW ISOFORM 25 Jul 2017
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData") ## NEW ISOFORM 25 Jul 2017

load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_paired_clinical.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
patients <- names(scores_list)
clin_data <- clin_data[rownames(clin_data) %in% patients,]


source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R")
source("~/Dropbox/Splice-n-of-1-pathways/Code/target_functions.R")

tmp_data <- scores_list[[1]]
(tmp_pat <- names(scores_list)[1])



############################################################
## 1. Panel A

\begin{table}[ht]
\centering
\begin{tabular}{lrrrrrl}
  \hline
Pathway description & Odds ratio & FDR & Critical OR & Alt splice \\ 
  \hline
Asthma & 4.05 & 0.00 & 1.93 & Yes \\ 
  Allograft rejection & 2.90 & 0.00 & 1.93 & Yes \\ 
  Taste transduction & 2.77 & 0.00 & 1.93 & Yes \\ 
  \vdots & \vdots & \vdots & \vdots & \vdots \\ 
  Wnt signaling pathway & 1.21 & 1.00 & 1.93 & No \\ 
  \vdots & \vdots & \vdots & \vdots & \vdots \\ 
   \hline
\end{tabular}
\end{table}

############################################################
## 2. Panel B

effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
rm_effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = T)
str(rm_effect_mat)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)

## filter to only those pathways found dysregulated in at least one patient
fil_effect_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "pathway_score")
fil_fdr_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "fdr_value")
str(fil_effect_mat) ## 101 pathways now
str(fil_fdr_mat)

b_data <- data.frame(effect_mat[1:6,1:6])

x <- xtable(b_data)
align(x) <- xalign(x)
display(x) <- xdisplay(x)
print.xtable(x, include.rownames = T)

% latex table generated in R 3.3.0 by xtable 1.8-2 package
% Mon Jul 31 13:15:33 2017
\begin{table}[ht]
\centering
\begin{tabular}{lrrrrrr}
  \hline
 Pat ID & Path 1 & Path 2 & Path 3 & Path 4 & \ldots & Path P \\ 
  \hline
22-4593 & 1.21 & 1.28 & 0.27 & 0.97 & \ldots & 0.92 \\ 
  22-4609 & 0.87 & 1.10 & 0.62 & 0.78 & \ldots & 0.29 \\ 
  22-5471 & 1.02 & 0.72 & 0.40 & 0.18 & \ldots& 0.43 \\ 
22-5472 & 1.03 & 0.88 & 0.00 & 0.98 & \ldots& 0.79 \\
  \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\ 
  22-5478 & 1.15 & 1.38 & 0.51 & 0.97 & \ldots& 0.55 \\ 
  22-5481 & 1.37 & 0.61 & 0.00 & 0.58 & \ldots& 0.46 \\ 
   \hline
\end{tabular}
\end{table}


############################################################
## 3. Panel C


% latex table generated in R 3.3.0 by xtable 1.8-2 package
% Mon Jul 31 13:15:33 2017
\begin{table}[ht]
\centering
\begin{tabular}{lrrrrrr}
  \hline
 Pat ID & Path 3 & Path 10 & Path 47 & Path 100 & \ldots & Path (P - F) \\ 
  \hline
22-4593 & 3.21 & 3.28 & 0.27 & 0.97 & \ldots & 2.92 \\ 
  22-4609 & 1.87 & 3.30 & 0.62 & 0.78 & \ldots & 0.29 \\ 
  22-5471 & 3.02 & 2.72 & 0.40 & 0.38 & \ldots& 0.43 \\ 
22-5472 & 3.03 & 0.88 & 0.00 & 0.98 & \ldots& 2.79 \\
  \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\ 
  22-5478 & 3.35 & 3.38 & 0.53 & 0.97 & \ldots& 3.55 \\ 
  22-5483 & 3.37 & 0.63 & 4.00 & 0.58 & \ldots& 0.46 \\ 
   \hline
\end{tabular}
\end{table}


############################################################
## 4. Panel D

(cancer_pathways <- get_cancer_pathways(scores_list))
(cancer_ids <- names(cancer_pathways))
(target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
(target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T))
target_data$TCGA <- "LUSC"
## add cluster to show patient heterogeneity
sort(target_data$OR)
(target_data$cluster <- cluster_pat(target_data$OR))

(panel_d <- ggplot(data = target_data, aes(x = factor(cluster), y = OR, color = factor(cluster))) +
     geom_boxplot() +
     theme(legend.position = "none") +
     background_grid(major = "xy", minor = "none") +
     ## ylim(c(0, 4)) +
     geom_jitter(width = 0.1, alpha = 0.5) +
     geom_hline(yintercept = 1) +
     xlab( "Cluster via KEGG pathway" ) +
     ylab( "Odds ratio of enriched alt splice")
     ##scale_color_discrete("TCGA", labels = c("BLCA (19)", "LUAD (58)", "LUSC (51)", "PRAD (52)", "THCA (59)", "UCEC (7)"))
     )

############################################################
## 5. Panel E, random survival curves

set.seed(444)
tmp_or <- fil_effect_mat[,sample(1:ncol(fil_effect_mat),1)]
(top_clust <- cluster_pat(tmp_or))

(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(panel_e <- top_surv$plot + labs(title = paste("Cluster via pathway X")))

############################################################
## 6. Panel F random pvalues

## mock pvalue
tmp_p <- 0.125
p_data <- data.frame(rand_pvalue = rand_pvalue)

(panel_f <- ggplot(data = p_data, aes(x = rand_pvalue)) +
     geom_histogram(bins = 120) +
     ## theme(legend.position = "none") +
     background_grid(major = "xy", minor = "none") +
     ## ylim(c(0, 4)) +
     ## geom_jitter(width = 0.1, alpha = 0.5) +
     geom_vline(xintercept = tmp_p, color = "red") +
     xlab( "Random p-values" ) +
     ylab( "Count")
     ##scale_color_discrete("TCGA", labels = c("BLCA (19)", "LUAD (58)", "LUSC (51)", "PRAD (52)", "THCA (59)", "UCEC (7)"))
     )

############################################################
## 7. Panel G, relevant pathways


\begin{table}[ht]
\centering
\begin{tabular}{lrrrrrl}
  \hline
Pathway description & Obs p-value & Emprical & FDR \\
\hline
Bladder cancer	& 0.00015 &	0.0005 &	0.0505 \\
Melanoma &	0.00240 &	0.0030 &	0.1178 \\
Staph &	0.00281 &	0.0030 &	0.1178 \\
Lung cancer & 0.00514	& 0.0070 &	0.1717 \\
   \hline
\end{tabular}
\end{table}


############################################################
## 8. Panel H

(cancer_pathways <- get_cancer_pathways(scores_list))
(cancer_ids <- names(cancer_pathways))
(target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
top_desc <- cancer_pathways[target]
(target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T))
target_data$TCGA <- "LUSC"
## add cluster to show patient heterogeneity
sort(target_data$OR)
(top_clust <- cluster_pat(target_data$OR))
names(top_clust) <- rownames(clin_data)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(panel_h <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))

############################################################
## 9. Panel I

############################################################
## 10. Combine and save

(fig2 <- plot_grid(NULL, NULL, NULL,
                   panel_d, panel_e, panel_f,
                   NULL, panel_h, NULL,
                   labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), ncol = 3))

save_plot("Figure2.pdf", fig2,
          ncol = 3, 
          nrow = 3,
          ## base_width = 5,
          base_aspect_ratio = 1.3
          )
