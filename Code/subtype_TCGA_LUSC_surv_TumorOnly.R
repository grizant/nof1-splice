## LUSC subtyping via using isoform RNAseq
## AG Schissler
## Created 29 Jul 2017

############################################################
## i. Load clinical data and results

load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_paired_clinical.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/lusc_iso_kegg_data.RData") ## NEW ISOFORM 25 Jul 2017

## load survival analysis functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R")

## explore clinical
table(clin_data$vitalstatus) ## 52 in clinical

## rename to make generic
iso_data <- lusc_iso_kegg_data
rm(lusc_iso_kegg_data)

## retrieve gene set annotations
kegg <- nof1::read_gene_set(file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt")
kegg_desc <- read.delim("~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", stringsAsFactors = F)
rownames(kegg_desc) <- as.character(kegg_desc$path_id)
str(kegg)
kegg_list <- split(kegg$symbol, kegg$path_id)

## filter to gene set size
size <- unlist(lapply(kegg_list, function(tmp_genes) {length(unique(tmp_genes))}))
kegg_list <- kegg_list[which(size >= 15 & size <= 500)]
(max_genes <- max(size <- unlist(lapply(kegg_list, function(tmp_genes) {length(unique(tmp_genes))}))))
(max_kegg <- which.max(size))

############################################################
## 1. Cluster on isoform expr and get surv empirical for each kegg pathway 

## use only tumor expression

iso_data <- iso_data[ , c(1,2, grep("-T", names(iso_data)))]
names(iso_data) <- gsub("-T$", "", names(iso_data))
## pathway_genes <- kegg_list[[1]]

## find patients and subset iso data
iso_patients <- names(iso_data[, -(1:2)])
length(iso_patients) ## 51 in iso data
clin_patients <- rownames(clin_data)
length(clin_patients)
## intersection
length(unique(iso_patients[iso_patients %in% clin_patients]))
iso_data <-  iso_data[ , c(1:2, which(names(iso_data) %in% clin_patients))]
ncol(iso_data) 

## test on the first and largest
system.time(print(get_empirical_pvalue(pathway_genes = kegg_list[[1]], iso_data, clin_data))) ## 13 sec for 65 genes
system.time(print(get_empirical_pvalue(pathway_genes = kegg_list[[max_kegg]], iso_data, clin_data))) ## 26 sec for 388 genes
 
## parallelize
num_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(num_cores, type = "FORK")

## export functions
parallel::clusterEvalQ(cl = cl, expr = source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R"))

## test on a subset
system.time(emp_pvalue <- unlist(parallel::parLapply(cl, X = kegg_list[sample(1:length(kegg_list), size = num_cores)], fun = get_empirical_pvalue, iso_data, clin_data)))
emp_pvalue

## run on all pathways
set.seed(44)
system.time(emp_pvalue <- unlist(parallel::parLapply(cl, X = kegg_list, fun = get_empirical_pvalue, iso_data, clin_data)))

str(emp_pvalue)

## stop Cluster
parallel::stopCluster(cl)

save(emp_pvalue, file = "~/Dropbox/Splice-n-of-1-pathways/Data/tumor_iso_surv_empirical_pvalues.RData")

############################################################
## 2. Explore the Emprical results

hist(emp_pvalue, breaks = 15)
emp_pvalue <- sort(emp_pvalue)
head(emp_pvalue)

## where is the target
target <- kegg_desc$path_id[grep("Non-small cell lung cancer", kegg_desc$description)]
emp_pvalue[target]
which(names(emp_pvalue) == target) ## 56

bh_pvalue <- p.adjust(emp_pvalue, "BH") ## nothing at FDR 50%
head(bh_pvalue, 56)
bh_pvalue[target]

threshold <- 0.05
any(emp_pvalue < threshold)
(hits <- names(emp_pvalue)[emp_pvalue < threshold])
sum(emp_pvalue < threshold)/length(emp_pvalue)

## retrive the hit pathway description
kegg_desc[hits,]

## try local FDR more direct comparison
transformed_pvalue <- qnorm(emp_pvalue)
eps <- 0.001
transformed_pvalue[is.infinite(transformed_pvalue)] <- max(transformed_pvalue[!is.infinite(transformed_pvalue)], na.rm = T) + eps

tmp_locfdr <- locfdr::locfdr(zz = transformed_pvalue, bre = ceiling(length(emp_pvalue)/8), df = 4, pct = 0, pct0 = 1/64, nulltype = 1, plot = 1)
## tmp_locfdr

## 0 hits by either 

######################################################
## 2.1 explore top 4 hits
library(ggplot2)
library(cowplot)

kegg_desc[hits[1:4], 2]

### hit 1
i <- 1
tmp_genes <- kegg_list[[hits[i]]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
one_clust <- top_clust
(top_desc <- kegg_desc[hits[i], 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s1 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))
    
### hit 2
i <- 2
tmp_genes <- kegg_list[[hits[i]]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
two_clust <- top_clust
(top_desc <- kegg_desc[hits[i], 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s2 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))

### hit 3
i <- 3
tmp_genes <- kegg_list[[hits[i]]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
three_clust <- top_clust
(top_desc <- kegg_desc[hits[i], 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s3 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))

### hit 4
i <- 4
tmp_genes <- kegg_list[[hits[i]]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
four_clust <- top_clust
(top_desc <- kegg_desc[hits[i], 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s4 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))

### target
tmp_genes <- kegg_list[[target]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
target_clust <- top_clust
(top_desc <- kegg_desc[target, 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s5 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))

##################
### combine together

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")

## Figure 4
surv_plots <- plot_grid(s1, s2, s3, s4,
                         labels=c("A", "B", "C", "D"), ncol = 2, align = "h")
save_plot("TumorIso_Figure4.pdf", surv_plots,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.1,
          )

save_plot("TumorIso_Target_pathway_surv.pdf", s5,)


############################################################
## 3. Compute observed p-values for local FDR

## parallelize
num_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(num_cores, type = "FORK")

## export functions
parallel::clusterEvalQ(cl = cl, expr = source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R"))

## test on a subset
system.time(obs_pvalues <- unlist(parallel::parLapply(cl, X = kegg_list[sample(1:length(kegg_list), size = num_cores)], fun = clust_surv_pvalue, iso_data, clin_data)))
obs_pvalue

## run on all pathways
set.seed(44)
system.time(obs_pvalues <- unlist(parallel::parLapply(cl, X = kegg_list, fun = clust_surv_pvalue, iso_data, clin_data)))

str(obs_pvalues)

## stop Cluster
parallel::stopCluster(cl)

save(obs_pvalues, file = "~/Dropbox/Splice-n-of-1-pathways/Data/tumor_iso_surv_obs_pvalues.RData")
 
############################################################
## 4. Explore the observed results

hist(obs_pvalues, breaks = 15)
obs_pvalues <- sort(obs_pvalues)
head(obs_pvalues)

bh_pvalue <- p.adjust(obs_pvalues, "BH") 
head(bh_pvalue, 4) ## nothing at FDR 20%

threshold <- 0.05
any(obs_pvalues < threshold)
(hits <- names(obs_pvalues)[obs_pvalues < threshold])
sum(obs_pvalues < threshold)/length(obs_pvalues)

## retrive the hit pathway description
kegg_desc[hits,]

## try local FDR more direct comparison
transformed_pvalue <- qnorm(obs_pvalues)
eps <- 0.001
transformed_pvalue[is.infinite(transformed_pvalue)] <- max(transformed_pvalue[!is.infinite(transformed_pvalue)], na.rm = T) + eps

tmp_locfdr <- locfdr::locfdr(zz = transformed_pvalue, bre = ceiling(length(obs_pvalues)/8), df = 4, pct = 0, pct0 = 1/64, nulltype = 1, plot = 1)
## tmp_locfdr

tmp_locfdr$fdr[which(names(obs_pvalues) == target)]

## 0 hits by either 

######################################################
## 2.1 explore top 4 hits
library(ggplot2)
library(cowplot)

kegg_desc[hits[1:4], 2]

### hit 1
i <- 1
tmp_genes <- kegg_list[[hits[i]]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
one_clust <- top_clust
(top_desc <- kegg_desc[hits[i], 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s1 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))
    
### hit 2
i <- 2
tmp_genes <- kegg_list[[hits[i]]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
two_clust <- top_clust
(top_desc <- kegg_desc[hits[i], 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s2 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))

### hit 3
i <- 3
tmp_genes <- kegg_list[[hits[i]]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
three_clust <- top_clust
(top_desc <- kegg_desc[hits[i], 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s3 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))

### hit 3
i <- 4
tmp_genes <- kegg_list[[hits[i]]]
tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
(top_clust <- cluster_pat(tmp_mat))
(obs_pvalue <- clust_surv_pvalue(tmp_genes = tmp_genes, iso_data, clin_data))
table(top_clust)
four_clust <- top_clust
(top_desc <- kegg_desc[hits[i], 2])
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T))
(s4 <- top_surv$plot + labs(title = paste("Cluster via '", top_desc, "'", sep="")))


##################
### combine together

setwd("~/Dropbox/splice-n-of-1-pathways/Figures")

## Figure 4
surv_plots <- plot_grid(s1, s2, s3, s4,
                         labels=c("A", "B", "C", "D"), ncol = 2, align = "h")
save_plot("TumorIso_locFDR_Figure4.pdf", surv_plots,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.1,
          )
