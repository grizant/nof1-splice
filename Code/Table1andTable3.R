## Validate TCGA all patients through aggregation of target pathway information
## AG Schissler
## Created 25 Jul 2017

##############################################################################
#### Source validation functions

source("~/Dropbox/Splice-n-of-1-pathways/Code/target_functions.R")

##load and aggregate manually
all_data <- NULL

## 1. BLCA
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(cancer_pathways <- get_cancer_pathways(scores_list))
(cancer_ids <- names(cancer_pathways))
(target <- names(cancer_pathways)[grep("Bladder", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
(to_return <- summarize_target_data(target_data))
rownames(to_return) <- "BLCA"
(all_data <- rbind(all_data, to_return))

## 2. THCA
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_THCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Thyroid", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
(to_return <- summarize_target_data(target_data))
rownames(to_return) <- "THCA"
(all_data <- rbind(all_data, to_return))

## 3. UCEC
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_UCEC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Endo", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
(to_return <- summarize_target_data(target_data))
rownames(to_return) <- "UCEC"
(all_data <- rbind(all_data, to_return))

## 4. PRAD
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_PRAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Prostate", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
(to_return <- summarize_target_data(target_data))
rownames(to_return) <- "PRAD"
(all_data <- rbind(all_data, to_return))

## 5. LUSC
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
(to_return <- summarize_target_data(target_data))
rownames(to_return) <- "LUSC"
(all_data <- rbind(all_data, to_return))

## ## 5.2 LUSC EEv2
## load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_EEv2_KEGG_29july2017.RData")
## (target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
## target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## (to_return <- summarize_target_data(target_data))
## rownames(to_return) <- "LUSCv2"
## (all_data <- rbind(all_data, to_return))
## 

## 6. LUAD
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
(target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
(to_return <- summarize_target_data(target_data))
rownames(to_return) <- "LUAD"
(all_data <- rbind(all_data, to_return))

## Table 3
all_data[order(all_data$target_capture_rate),]

## ## READ (drop as this target is not a perfect match)
## load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_READ_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## (target <- names(cancer_pathways)[grep("Colorectal", cancer_pathways)])
## target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## (to_return <- summarize_target_data(target_data))
## rownames(to_return) <- "READ"
## (all_data <- rbind(all_data, to_return))
## 
## ## COAD
## load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_COAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## (target <- names(cancer_pathways)[grep("Colorectal", cancer_pathways)])
## target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## (to_return <- summarize_target_data(target_data))
## rownames(to_return) <- "COAD"
## (all_data <- rbind(all_data, to_return))
## 

##### Table 1

## get clinical data

## LUSC
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## LUAD
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUAD_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## PRAD
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_PRAD_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## THCA
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_THCA_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## UCEC
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_UCEC_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## BLCA
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

